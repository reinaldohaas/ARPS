!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADVPTFCT                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advptfct(nx,ny,nz,dtbig1,ptprt,u,v,w,wcont,                  &
           rhostr,rhostri,ptbar,mapfct,j3,j3inv,                        &
           ptadv,                                                       &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                     &
           ustr,vstr,wstr,tem1_0,tem2_0,tem3_0,mp_tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the advection of total potential temperature
!  using flux-correctd transport scheme.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue and Yifeng Tang
!  1/6/93.
!
!  10/12/1996 (M. Xue)
!  Code clean up. prprt*rhostr is now passed into advctsfct.
!
!  7/10/1997 (Fanyou Kong -- CMRP)
!  Include MPDCD simple positive definite advection scheme
!  (sadvopt = 5) through subroutine 'advmpdcd'
!
!  11/20/1997 (Fanyou Kong -- CMRP)
!  Added map factor into both FCT and MPDCD schemes
!
!  7/17/1998 (Ming Xue)
!  Significant modifications to call a new version of advctsfct.
!
!  10/5/1998 (Dan Weber)
!  Added rhostri, mapfct 7-8 array's.
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
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature at all time levels (K)
!
!    rhostr   Base state air density times j3 (kg/m**3)
!    rhostri  Inverse base state air density times j3 (kg/m**3)
!    ptbar    Base state potential temperature (K)
!    mapfct   Map factor for scalar point
!
!  OUTPUT:
!
!    ptadv    Advection term of potential temperature eqn (kg/m**3)*K/s
!
!  WORK ARRAYS:
!
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!    ustr     Work array
!    vstr     Work array
!    wstr     Work array
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
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INCLUDE 'timelvls.inc'

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL :: dtbig1

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: u     (nx,ny,nz,nt)  ! u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant w-velocity (m/s)
!
  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3 (kg/m**3)
  REAL :: rhostri(nx,ny,nz)    ! Inv. base state air density times j3 (kg/m**3)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: mapfct(nx,ny,8)      ! Map factor for scalar point
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ).
  REAL :: j3inv (nx,ny,nz)     ! 1/J3.
!
  REAL :: ptadv(nx,ny,nz)      ! Advection of potential temperature
!
!
  REAL :: tem1  (nx,ny,nz)     ! temporary work array
  REAL :: tem2  (nx,ny,nz)     ! temporary work array
  REAL :: tem3  (nx,ny,nz)     ! temporary work array
  REAL :: tem4  (nx,ny,nz)     ! temporary work array
  REAL :: tem5  (nx,ny,nz)     ! temporary work array
  REAL :: tem6  (nx,ny,nz)     ! temporary work array
  REAL :: tem7  (nx,ny,nz)     ! temporary work array
  REAL :: tem8  (nx,ny,nz)     ! temporary work array
  REAL :: ustr  (nx,ny,nz)     ! Work array holding u*rhostr/mapfct
  REAL :: vstr  (nx,ny,nz)     ! Work array holding v*rhostr/mapfct
  REAL :: wstr  (nx,ny,nz)     ! Work array holding wcont*rhostr

  REAL :: tem1_0(0:nx,0:ny,0:nz) ! Work array
  REAL :: tem2_0(0:nx,0:ny,0:nz) ! Work array
  REAL :: tem3_0(0:nx,0:ny,0:nz) ! Work array

  REAL :: mp_tem(MAX(nx+1,ny+1)*(nz+1))  ! Temporary message passing array.

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: ptmin, ptmin1, ptmax1, ptmin2, ptmax2, tema

  INTEGER :: advdiv     ! = 0, neglect anelastic correction in FCT advection
                        ! = 1, include anelastic correction in FCT advection
!  INTEGER :: mptag1,mptag2
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'       ! Global model control parameters
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'bndry.inc'          ! Boundary Condition Flags
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  advdiv    = 1  ! Include the anelastic correction term? (recommended!)

  IF( fctadvptprt == 1 .AND. sadvopt == 5 ) fctadvptprt = 2

  IF( ptsmlstp == 1 ) fctadvptprt = 1  ! Do not reset this one

  IF( fctadvptprt == 0 ) THEN

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1_0(i,j,k)=ptprt(i,j,k,1)+ptbar(i,j,k)
          tem2_0(i,j,k)=ptprt(i,j,k,2)+ptbar(i,j,k)
        END DO
      END DO
    END DO

  ELSE IF( fctadvptprt == 1 ) THEN

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1_0(i,j,k)=ptprt(i,j,k,1)
          tem2_0(i,j,k)=ptprt(i,j,k,2)
        END DO
      END DO
    END DO

  ELSE IF( fctadvptprt == 2 ) THEN

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=ptprt(i,j,k,1)+ptbar(i,j,k)
          tem2(i,j,k)=ptprt(i,j,k,2)+ptbar(i,j,k)
        END DO
      END DO
    END DO

    CALL a3dmax0(tem1,1,nx,1,nx-1,1,ny,1,ny-1,                          &
                 1,nz,1,nz-1,ptmax1,ptmin1)
    CALL a3dmax0(tem2,1,nx,1,nx-1,1,ny,1,ny-1,                          &
                 1,nz,1,nz-1,ptmax2,ptmin2)
    ptmin = AMIN1(ptmin1,ptmin2)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1_0(i,j,k)=tem1(i,j,k)-ptmin
          tem2_0(i,j,k)=tem2(i,j,k)-ptmin
        END DO
      END DO
    END DO

  END IF

!  IF (mp_opt > 0) THEN
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendextew(tem1_0,nx,ny,nz,ebc,wbc,mptag1,mp_tem)
!    CALL mpsendextew(tem2_0,nx,ny,nz,ebc,wbc,mptag2,mp_tem)
!    CALL mprecvextew(tem1_0,nx,ny,nz,ebc,wbc,mptag1,mp_tem)
!    CALL mprecvextew(tem2_0,nx,ny,nz,ebc,wbc,mptag2,mp_tem)
!    CALL mpsendextns(tem1_0,nx,ny,nz,nbc,sbc,mptag1,mp_tem)
!    CALL mpsendextns(tem2_0,nx,ny,nz,nbc,sbc,mptag2,mp_tem)
!    CALL mprecvextns(tem1_0,nx,ny,nz,nbc,sbc,mptag1,mp_tem)
!    CALL mprecvextns(tem2_0,nx,ny,nz,nbc,sbc,mptag2,mp_tem)
!  END IF
  CALL acct_interrupt(bc_acct)
  CALL extndsbc(tem1_0,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)
  CALL extndsbc(tem2_0,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)
  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
!  Calculate ustr, vstr and wstr.
!
!-----------------------------------------------------------------------
!

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem3_0(i,j,k)=rhostr(i,j,k)
      END DO
    END DO
  END DO

!  IF (mp_opt > 0) THEN
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendextew(tem3_0,nx,ny,nz,ebc,wbc,mptag1,mp_tem)
!    CALL mprecvextew(tem3_0,nx,ny,nz,ebc,wbc,mptag1,mp_tem)
!    CALL mpsendextns(tem3_0,nx,ny,nz,nbc,sbc,mptag1,mp_tem)
!    CALL mprecvextns(tem3_0,nx,ny,nz,nbc,sbc,mptag1,mp_tem)
!  END IF
  CALL acct_interrupt(bc_acct)
  CALL extndsbc(tem3_0,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)
  CALL acct_stop_inter

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        ustr(i,j,k)=u(i,j,k,2)*(tem3_0(i-1,j,k)+tem3_0(i,j,k))          &
                   *mapfct(i,j,5)*0.5
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        vstr(i,j,k)=v(i,j,k,2)*(tem3_0(i,j-1,k)+tem3_0(i,j,k))          &
                   *mapfct(i,j,6)*0.5
      END DO
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx-1
        wstr(i,j,k)=wcont(i,j,k)                                        &
                   *(tem3_0(i,j,k-1)+tem3_0(i,j,k))*0.5
      END DO
    END DO
  END DO

  IF(advdiv == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Add the anelastic divergence correction term:
! -[m**2{d(rhostr*u/m)/dx+ d(rhostr*v/m)/dy}+ d(rhostr*wcont)/dz]*ptprt.
!
!-----------------------------------------------------------------------
!
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          ptadv(i,j,k)= - tem2_0(i,j,k)*(                               &
               ((ustr(i+1,j,k)-ustr(i,j,k))*dxinv                       &
               +(vstr(i,j+1,k)-vstr(i,j,k))*dyinv)*mapfct(i,j,7)        &
               +(wstr(i,j,k+1)-wstr(i,j,k))*dzinv )
        END DO
      END DO
    END DO

  ELSE

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          ptadv(i,j,k)= 0.0
        END DO
      END DO
    END DO

  END IF

  IF (sadvopt == 4) THEN

    IF(fctadvptprt == 1.AND.ptsmlstp /= 1 ) THEN  ! Advection of ptbar

      tema = 0.25/dz

      DO k=2,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem2(i,j,k)=(rhostr(i,j,k  )*j3inv(i,j,k  )                 &
                        +rhostr(i,j,k-1)*j3inv(i,j,k-1))*w(i,j,k,2)     &
                       *(ptbar(i,j,k)-ptbar(i,j,k-1))*tema
          END DO
        END DO
      END DO

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            ptadv(i,j,k)=ptadv(i,j,k)+tem2(i,j,k)+tem2(i,j,k+1)
          END DO
        END DO
      END DO

    END IF

    CALL advctsfct(nx,ny,nz,nt,dtbig1,ustr,vstr,wstr,                   &
         mapfct(:,:,1),tem1_0,tem2_0,tem3_0,rhostr,rhostri,             &
         ptadv,                                                         &
         tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)

  ELSE IF(sadvopt == 5) THEN

    IF(fctadvptprt /= 0 ) THEN
      CALL arpsstop(                                                    &
          'Subourinte ADVMPDCD should never be used to advect '//       &
          'non-positive definite variable PTPRT.'//                     &
          'Job Stopped in subroutine ADVPTFCT.',1)
    END IF

    CALL advmpdcd(nx,ny,nz,nt,dtbig1,                                   &
         ustr,vstr,wstr,mapfct(:,:,1),tem1_0,tem2_0, rhostr,            &
         ptadv,                                                         &
         tem1,tem2,tem3,tem4, tem3_0, mp_tem)

  END IF

  RETURN
END SUBROUTINE advptfct

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADVQFCT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advqfct(nx,ny,nz,dtbig1,q,u,v,wcont, ustr,vstr,wstr,         &
           rhostr,rhostri,mapfct,j3,qadv,                               &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                     &
           tem1_0,tem2_0,tem3_0)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the advection of water substance
!  using flux-correctd transport scheme.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yifeng Tang and M. Xue
!  1/6/93.
!
!  10/12/1996 (M. Xue)
!  Code clean up. q*rhostr is now passed into advctsfct.
!
!  7/10/1997 (Fanyou Kong -- CMRP)
!  Include MPDCD simple positive definite advection scheme
!  (sadvopt = 5) through subroutine 'advmpdcd'
!
!  11/20/1997 (Fanyou Kong -- CMRP)
!  Added map factor into FCT and MPDCD schemes
!
!  7/17/1998 (Ming Xue)
!  Significant modifications to call a new version of advctsfct.
!
!  9/28/1998 (Dan Weber)
!  Moved calculation of u,v,wstr out of this subroutine and added
!  rhostri and mapfct 7-8.
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
!    wcont    Contravariant vertical velocity (m/s)
!    ustr     Appropriate ustr for advection
!    vstr     Appropriate ustr for advection
!    wstr     Appropriate ustr for advection
!    q        Water substance at all time levels (K)
!
!    rhostr   Base state air density times j3 (kg/m**3)
!    rhostri  Inverse base state air density times j3 (kg/m**3)
!    mapfct   Map factor for scalar point
!
!  OUTPUT:
!
!    qadv    Advection term of water substance eqn (kg/m**3)*(kg/kg)/s
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
!    tem1_0   Temporary work array.
!    tem2_0   Temporary work array.
!    tem3_0   Temporary work array.
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
  INCLUDE 'timelvls.inc'
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL :: dtbig1

  REAL :: q     (nx,ny,nz,nt)  ! Water substance  (kg/kg)
  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ustr  (nx,ny,nz)     ! u*rhostr/mapfct
  REAL :: vstr  (nx,ny,nz)     ! v*rhostr/mapfct
  REAL :: wstr  (nx,ny,nz)     ! wcont*rhostr
!
  REAL :: rhostr(nx,ny,nz)     ! Base state air density times j3 (kg/m**3)
  REAL :: rhostri(nx,ny,nz)    ! Inv. base state air density times j3 (kg/m**3)
  REAL :: mapfct(nx,ny,8)      ! Map factor for scalar point
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ).
!
  REAL :: qadv  (nx,ny,nz)     ! Advection of water substance
!
  REAL :: tem1  (nx,ny,nz)     ! temporary work array
  REAL :: tem2  (nx,ny,nz)     ! temporary work array
  REAL :: tem3  (nx,ny,nz)     ! temporary work array
  REAL :: tem4  (nx,ny,nz)     ! temporary work array
  REAL :: tem5  (nx,ny,nz)     ! temporary work array
  REAL :: tem6  (nx,ny,nz)     ! temporary work array
  REAL :: tem7  (nx,ny,nz)     ! temporary work array
  REAL :: tem8  (nx,ny,nz)     ! temporary work array
  REAL :: tem1_0(0:nx,0:ny,0:nz) ! automatic work array
  REAL :: tem2_0(0:nx,0:ny,0:nz) ! automatic work array
  REAL :: tem3_0(0:nx,0:ny,0:nz) ! automatic work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: advdiv

!  INTEGER :: mptag1,mptag2
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'       ! Global model control parameters
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'bndry.inc'         ! Boundary Condition Flags
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  advdiv = 1  ! Include anelastic correction term in advection
!
!-----------------------------------------------------------------------
!
!  Calculate the q(tfurture) with FCT correction.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1_0(i,j,k)=q(i,j,k,1)
        tem2_0(i,j,k)=q(i,j,k,2)
      END DO
    END DO
  END DO

!  IF (mp_opt > 0) THEN
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendextew(tem1_0,nx,ny,nz,ebc,wbc,mptag1,tem1)
!    CALL mpsendextew(tem2_0,nx,ny,nz,ebc,wbc,mptag2,tem1)
!    CALL mprecvextew(tem1_0,nx,ny,nz,ebc,wbc,mptag1,tem1)
!    CALL mprecvextew(tem2_0,nx,ny,nz,ebc,wbc,mptag2,tem1)
!    CALL mpsendextns(tem1_0,nx,ny,nz,nbc,sbc,mptag1,tem1)
!    CALL mpsendextns(tem2_0,nx,ny,nz,nbc,sbc,mptag2,tem1)
!    CALL mprecvextns(tem1_0,nx,ny,nz,nbc,sbc,mptag1,tem1)
!    CALL mprecvextns(tem2_0,nx,ny,nz,nbc,sbc,mptag2,tem1)
!  END IF
  CALL acct_interrupt(bc_acct)
  CALL extndsbc(tem1_0,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)
  CALL extndsbc(tem2_0,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)
  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
!  Calculate ustr, vstr and wstr.
!
!-----------------------------------------------------------------------
!

  IF( advdiv == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Add the anelastic divergence correction term:
!  -[d(rhostr*u)/dx+ d(rhostr*v)/dy+ d(rhostr*wcont)/dz]*q
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          qadv(i,j,k)= - tem2_0(i,j,k)*(                                &
              ((ustr(i+1,j,k)-ustr(i,j,k))*dxinv                        &
              +(vstr(i,j+1,k)-vstr(i,j,k))*dyinv)*mapfct(i,j,7)         &
              +(wstr(i,j,k+1)-wstr(i,j,k))*dzinv)
        END DO
      END DO
    END DO

  ELSE

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          qadv(i,j,k)= 0.0
        END DO
      END DO
    END DO

  END IF

  IF(sadvopt == 4) THEN

    CALL advctsfct(nx,ny,nz,nt,dtbig1,ustr,vstr,wstr,    &
         mapfct(1,1,7),tem1_0,tem2_0,tem3_0,rhostr,rhostri,             &
         qadv,                                                          &
         tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)

  ELSE IF(sadvopt == 5) THEN  ! MPDCD positive definite advection

    CALL advmpdcd(nx,ny,nz,nt,dtbig1,ustr,vstr,wstr,           &
         mapfct(1,1,7),tem1_0,tem2_0,rhostr,                &
         qadv,                                                          &
         tem1,tem2,tem3,tem4,tem3_0,tem5)

  END IF

  RETURN
END SUBROUTINE advqfct
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADVCTSFCT                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advctsfct(nx,ny,nz,nt,dtbig1,                 &
           ustr,vstr,wstr,mapfct2, s1,s2, s3,rhostr,rhostri,            &
           sadv,                                                        &
           fluxx,fluxy,fluxz,cx,cy,cz,tem1,tem2)
!
!----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Integrate scalar field forward on time step using Zalesak's
!  multi-dimensional version of FCT.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue. Adaptied for ARPS by Y.F. Tang.
!  11/29/1992.
!
!  3/10/1993 (M. Xue)
!  Code cleanup.
!
!  11/20/1997 (Fanyou Kong -- CMRP)
!  Added map factor into FCT schemes
!
!  7/17/1998 (Ming Xue)
!
!  The formulation of advective fluxes are now consistent with
!  the regular advection.  Array sadv is now accumulative.
!
!  10/5/1998 (Dan Weber)
!  Added rhostri array.
!
!  2000/09/15 (Gene Bassett)
!  Removed dx,dy,dz,dtbig & fctorderopt from argument list.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nt       number of the time level
!    dtbig1   Time step size
!
!    ustr     u*rhostr
!    vstr     v*rhostr
!    wstr     wcont*rhostr
!    mapfct2  map factor for the scalar ** 2
!    s1       scalar field at a given time n-1 level (scalar units)
!    s2       scalar field at a given time n level (scalar units)
!    rhostr   rhobar*j3
!    rhostri  inverse rhobar*j3
!    sadv     Part of the advection term for scalar s
!
!  OUTPUT:
!
!    sadv     Advection for scalar s
!    s1,s2 and s3 modified.
!
!  WORK ARRAYS:
!
!    fluxx     Temporary work array.
!    fluxy     Temporary work array.
!    fluxz     Temporary work array.
!    cx        Temporary work array.
!    cy        Temporary work array.
!    cz        Temporary work array.
!    tem1      Temporary work array.
!    tem2      Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz
  INTEGER :: nt              ! The no. of t-levels of t-dependent arrays.
  REAL :: dtbig1

  REAL :: ustr  (nx,ny,nz)   ! u*rhostr
  REAL :: vstr  (nx,ny,nz)   ! v*rhostr
  REAL :: wstr  (nx,ny,nz)   ! wcont*rhostr
  REAL :: mapfct2(nx,ny)     ! map factor for scalar **2
  REAL :: s1    (0:nx,0:ny,0:nz)! Time lev. 1 of a scalar field to be advected.
  REAL :: s2    (0:nx,0:ny,0:nz)! Time lev. 2 of a scalar field to be advected.
  REAL :: s3    (0:nx,0:ny,0:nz)! Time lev. 3 of a scalar field to be advected.
  REAL :: rhostr(nx,ny,nz)   ! rhostr*j3
  REAL :: rhostri(nx,ny,nz)  ! Inverse rhostr*j3
  REAL :: sadv  (nx,ny,nz)   ! Advection of scalar s
!
  REAL :: fluxx (nx,ny,nz)   ! temporary work array
  REAL :: fluxy (nx,ny,nz)   ! temporary work array
  REAL :: fluxz (nx,ny,nz)   ! temporary work array
  REAL :: cx    (nx,ny,nz)   ! temporary work array
  REAL :: cy    (nx,ny,nz)   ! temporary work array
  REAL :: cz    (nx,ny,nz)   ! temporary work array
  REAL :: tem1  (nx,ny,nz)   ! temporary work array
  REAL :: tem2  (nx,ny,nz)   ! temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: tem, rdxyz,dxyz, tema,rdt 

!  INTEGER :: mptag1,mptag2
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'       ! Global model control parameters
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'bndry.inc'         ! Boundary Condition Flags
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  dxyz=dx*dy*dz
  rdxyz=1./(dxyz)

!-----------------------------------------------------------------------
!
!  Compute the second-order advective fluxes
!
!-----------------------------------------------------------------------
!
  CALL advflxs(nx,ny,nz,nt,dx,dy,dz,dtbig1,                             &
               ustr,vstr,wstr,s2,fctorderopt,                           &
               fluxx,fluxy,fluxz,                                       &
               s3,tem1)  ! Attention: s3 is used as a work array here!
!
!-----------------------------------------------------------------------
!
!  Compute the temperary scalar field at tfuture
!
!-----------------------------------------------------------------------
!
  tem = 2.0 * rdxyz
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        s3(i,j,k)=s1(i,j,k)-((fluxx(i+1,j,k)-fluxx(i,j,k)+              &
                              fluxy(i,j+1,k)-fluxy(i,j,k))*mapfct2(i,j)+ &
                              fluxz(i,j,k+1)-fluxz(i,j,k))*tem          &
                  *rhostri(i,j,k)
      END DO
    END DO
  END DO

!  IF (mp_opt > 0) THEN
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendextew(s3,nx,ny,nz,ebc,wbc,mptag1,tem1)
!    CALL mprecvextew(s3,nx,ny,nz,ebc,wbc,mptag1,tem1)
!    CALL mpsendextns(s3,nx,ny,nz,nbc,sbc,mptag1,tem1)
!    CALL mprecvextns(s3,nx,ny,nz,nbc,sbc,mptag1,tem1)
!  END IF
  CALL acct_interrupt(bc_acct)
  CALL extndsbc(s3,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)
  CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Trapezoidal step:
!
!-----------------------------------------------------------------------

  DO k=0,nz
    DO j=0,ny
      DO i=0,nx
        s3(i,j,k)=(s3(i,j,k)+s2(i,j,k))*0.5
      END DO
    END DO
  END DO

  CALL advflxs(nx,ny,nz,nt,dx,dy,dz,dtbig,                              &
               ustr,vstr,wstr,s3,fctorderopt,                           &
               fluxx,fluxy,fluxz,                                       &
               s1,tem1)  ! Attention: s1 is used as a work array here!
!
!-----------------------------------------------------------------------
!
!  Low order scheme: Donnor Cell.
!  Calculate the donnor-cell fluxes
!
!-----------------------------------------------------------------------

  tem = dtbig*dy*dz*0.5

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        cx(i,j,k)=tem*                                                  &
            ((ustr(i,j,k)+ABS(ustr(i,j,k)))*s2(i-1,j,k)                 &
            +(ustr(i,j,k)-ABS(ustr(i,j,k)))*s2(i  ,j,k))
        fluxx(i,j,k)=fluxx(i,j,k)-cx(i,j,k)
      END DO
    END DO
  END DO
!
  tem = dtbig*dx*dz*0.5

  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        cy(i,j,k)=tem*                                                  &
            ((vstr(i,j,k)+ABS(vstr(i,j,k)))*s2(i,j-1,k)                 &
            +(vstr(i,j,k)-ABS(vstr(i,j,k)))*s2(i,j  ,k))
        fluxy(i,j,k)=fluxy(i,j,k)-cy(i,j,k)
      END DO
    END DO
  END DO
!
  tem = dtbig*dx*dy*0.5

  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx-1
        cz(i,j,k)=tem*                                                  &
            ((wstr(i,j,k)+ABS(wstr(i,j,k)))*s2(i,j,k-1)                 &
            +(wstr(i,j,k)-ABS(wstr(i,j,k)))*s2(i,j,k  ))
        fluxz(i,j,k)=fluxz(i,j,k)-cz(i,j,k)
      END DO
    END DO
  END DO
!
!  Compute the low order time advanced solution for scalar
!
  rdt=1.0/dtbig

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tema = ((cx(i+1,j,k)-cx(i,j,k)                                  &
                +cy(i,j+1,k)-cy(i,j,k))*mapfct2(i,j)                    &
                +cz(i,j,k+1)-cz(i,j,k))*rdxyz
        s3(i,j,k)=s2(i,j,k)-tema*rhostri(i,j,k)

        sadv(i,j,k) = sadv(i,j,k)+ tema * rdt

      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Determine the flux limitor by calling fctlim.
!  fctlim returns the limited/corrected anti-diffusive fluxes
!  which are added to the lower order fluxes.
!
!c#######################################################################

!  IF (mp_opt > 0) THEN
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendextew(s2,nx,ny,nz,ebc,wbc,mptag1,tem1)
!    CALL mpsendextew(s3,nx,ny,nz,ebc,wbc,mptag2,tem1)
!    CALL mprecvextew(s2,nx,ny,nz,ebc,wbc,mptag1,tem1)
!    CALL mprecvextew(s3,nx,ny,nz,ebc,wbc,mptag2,tem1)
!    CALL mpsendextns(s2,nx,ny,nz,nbc,sbc,mptag1,tem1)
!    CALL mpsendextns(s3,nx,ny,nz,nbc,sbc,mptag2,tem1)
!    CALL mprecvextns(s2,nx,ny,nz,nbc,sbc,mptag1,tem1)
!    CALL mprecvextns(s3,nx,ny,nz,nbc,sbc,mptag2,tem1)
!  END IF
  CALL acct_interrupt(bc_acct)
  CALL extndsbc(s2,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)
  CALL extndsbc(s3,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)
  CALL acct_stop_inter

  CALL fctlim(nx,ny,nz,dxyz,s3,s2,rhostr,mapfct2,fluxx,fluxy,fluxz,     &
              cx,cy,cz,tem1,tem2)

!-----------------------------------------------------------------------
!
!  Add the flux-limitted adtidifusive fluxes.
!  Note the flux nomal to the lateral boundaries are not applied.
!
!c#######################################################################

  tem = 1.0/(dtbig*dx*dy*dz)

  DO k=2,nz-2

    DO j=1,ny-1
      DO i=2,nx-2
        sadv(i,j,k)=sadv(i,j,k)+(fluxx(i+1,j,k)-fluxx(i,j,k))           &
                   *mapfct2(i,j)*tem
      END DO
    END DO

    DO j=2,ny-2
      DO i=1,nx-1
        sadv(i,j,k)=sadv(i,j,k)+(fluxy(i,j+1,k)-fluxy(i,j,k))           &
                   *mapfct2(i,j)*tem
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        sadv(i,j,k)=sadv(i,j,k)+(fluxz(i,j,k+1)-fluxz(i,j,k))*tem
      END DO
    END DO

  END DO

  RETURN
END SUBROUTINE advctsfct
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADVFLXS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advflxs(nx,ny,nz,nt,dx,dy,dz,dt,                             &
           ustr,vstr,wstr,s,advodropt,                                  &
           fluxx,fluxy,fluxz,                                           &
           tem1_0,mp_tem)
!
!----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate 2nd- and 4th-order fluxes for the advection of a scalar.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue.
!  7/18/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nt       number of the time level
!    dx       Grid interval in x direction
!    dy       Grid interval in x direction
!    dz       Grid interval in x direction
!    dt       Time step size
!
!    ustr     u*rhostr
!    vstr     v*rhostr
!    wstr     wcont*rhostr
!    s        Scalar whose advective fluxes are to be found
!
!  OUTPUT:
!
!    fluxx    Advective flux in x direction.
!    fluxy    Advective flux in y direction.
!    fluxz    Advective flux in z direction.
!
!  Work array:
!    tem1_0   Work array
!    mp_tem   Work array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz
  INTEGER :: nt              ! The no. of t-levels of t-dependent arrays.
  REAL :: dx,dy,dz,dt

  REAL :: ustr  (nx,ny,nz)   ! u*rhostr
  REAL :: vstr  (nx,ny,nz)   ! v*rhostr
  REAL :: wstr  (nx,ny,nz)   ! wcont*rhostr
  REAL :: s(0:nx,0:ny,0:nz)  ! Scalar whose advective fluxes are to be found
!
  REAL :: fluxx (nx,ny,nz)   ! Advective flux in x direction.
  REAL :: fluxy (nx,ny,nz)   ! Advective flux in y direction.
  REAL :: fluxz (nx,ny,nz)   ! Advective flux in z direction.

  REAL :: tem1_0(0:nx,0:ny,0:nz) ! Work array
  REAL :: mp_tem(MAX(nx+1,ny+1)*(nz+1))  ! Temporary message passing array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!  REAL :: tem, rdxyz,dxyz, tem16, tem43
  REAL :: tem, tem16, tem43

  INTEGER :: advodropt

!  INTEGER :: mptag
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'     ! the boundary flags
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  dxyz=dx*dy*dz
!  rdxyz=1./(dxyz)

  tem43 = 4.0/3.0
  tem16 = 1.0/6.0
!
!-----------------------------------------------------------------------
!
!  Compute the second-order advective fluxes
!
!-----------------------------------------------------------------------
!
  tem = dt*dy*dz*0.5
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        fluxx(i,j,k)=tem*ustr(i,j,k)*(s(i,j,k)+s(i-1,j,k))
      END DO
    END DO
  END DO

  IF( advodropt == 2) THEN  ! For fourth order fluxes

    tem = dt*dy*dz*0.25
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1_0(i,j,k)=                                                &
              tem*(ustr(i,j,k)+ustr(i+1,j,k))*(s(i+1,j,k)+s(i-1,j,k))
        END DO
      END DO
    END DO

!    IF (mp_opt > 0) THEN
!      CALL acct_interrupt(mp_acct)
!      CALL mpsendextew(tem1_0,nx,ny,nz,ebc,wbc,mptag,mp_tem)
!      CALL mprecvextew(tem1_0,nx,ny,nz,ebc,wbc,mptag,mp_tem)
!    END IF
    CALL acct_interrupt(bc_acct)
    CALL extndsbc(tem1_0,nx,ny,nz,1,ebc,wbc,0,0,0,0)
    CALL acct_stop_inter

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          fluxx(i,j,k)=tem43*fluxx(i,j,k)                               &
                      -tem16*(tem1_0(i-1,j,k)+tem1_0(i,j,k))
        END DO
      END DO
    END DO

  END IF

  tem = dt*dx*dz*0.5
  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        fluxy(i,j,k)=tem*vstr(i,j,k)*(s(i,j,k)+s(i,j-1,k))
      END DO
    END DO
  END DO

  IF( advodropt == 2) THEN  ! For fourth order fluxes

    tem = dt*dx*dz*0.25
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1_0(i,j,k)=                                                &
              tem*(vstr(i,j,k)+vstr(i,j+1,k))*(s(i,j+1,k)+s(i,j-1,k))
        END DO
      END DO
    END DO

!    IF (mp_opt > 0) THEN
!      CALL acct_interrupt(mp_acct)
!      CALL mpsendextns(tem1_0,nx,ny,nz,nbc,sbc,mptag,mp_tem)
!      CALL mprecvextns(tem1_0,nx,ny,nz,nbc,sbc,mptag,mp_tem)
!    END IF
    CALL acct_interrupt(bc_acct)
    CALL extndsbc(tem1_0,nx,ny,nz,1,0,0,nbc,sbc,0,0)
    CALL acct_stop_inter

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1
          fluxy(i,j,k)=tem43*fluxy(i,j,k)                               &
                      -tem16*(tem1_0(i,j-1,k)+tem1_0(i,j,k))
        END DO
      END DO
    END DO

  END IF

  tem = dt*dx*dy*0.5
  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx-1
        fluxz(i,j,k)=tem*wstr(i,j,k)*(s(i,j,k)+s(i,j,k-1))
      END DO
    END DO
  END DO

  IF( advodropt == 2) THEN  ! For fourth order fluxes

    tem = dt*dx*dy*0.25
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1_0(i,j,k)=                                                &
              tem*(wstr(i,j,k)+wstr(i,j,k+1))*(s(i,j,k+1)+s(i,j,k-1))
        END DO
      END DO
    END DO

    CALL acct_interrupt(bc_acct)
    CALL extndsbc(tem1_0,nx,ny,nz,1,0,0,0,0,tbc,bbc)
    CALL acct_stop_inter

    DO k=1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          fluxz(i,j,k)=tem43*fluxz(i,j,k)                               &
                      -tem16*(tem1_0(i,j,k-1)+tem1_0(i,j,k))
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE advflxs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADVMPDCD                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Coastal Meteorology Research Program (CMRP)      ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advmpdcd(nx,ny,nz,nt,dtbig1,                                 &
           ustr,vstr,wstr,mapfct2,s1,s2,rhostr,                         &
           sadv,                                                        &
           fluxx,fluxy,fluxz,fout,bout,mp_tem)
!
!----------------------------------------------------------------------
!
!  PURPOSE:
!
!  MPDCD (Multidimensional Positive Definite Centered Difference
!  Scheme) positive definite advection for scalars
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Fanyou Kong (CMRP)
!  07/10/1997.
!
!  11/20/1997 (Fanyou Kong -- CMRP)
!  - Added map factor into MPDCD schemes
!  - Added Loop 50 to record possible negative scalar values in "s1"
!    and then force them positive to prevent possible floating
!    point error, but report the situation in output file
!
!  7/17/1998 (Ming Xue)
!    Modified the formulation of fluxes. sadv is now accumulative.
!
!  9/12/1999 (Pengfei Zhang & Ming Xue)
!    Added a constant "eps" in the denominator of an equation to
!    prevent overflow.
!
!  2000/09/15 (Gene Bassett)
!  Removed dx,dy,dz & fctorderopt from argument list.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nt       number of the time level
!    dx       Grid interval in x direction
!    dy       Grid interval in x direction
!    dz       Grid interval in x direction
!    dtbig1   Time step size
!
!    ustr     u*rhostr
!    vstr     v*rhostr
!    wstr     wcont*rhostr
!
!    mapfct2  map factor for scalar **2
!    s1       scalar field at a given time n-1 level (scalar units)
!    s2       scalar field at a given time n level (scalar units)
!    rhostr   rhobar*j3
!    fctorderopt Option parameter for 2nd or 4th order high-order fluxes
!
!  OUTPUT:
!
!    sadv     advection term              [(kg/m**3) scalar unit/s]
!
!  WORK ARRAYS:
!
!    fluxx     Temporary work array.
!    fluxy     Temporary work array.
!    fluxz     Temporary work array.
!    bout      Temporary work array.
!    fout      Temporary work array.
!    mp_tem    Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz
  INTEGER :: nt              ! The no. of t-levels of t-dependent arrays.
  REAL :: dtbig1

  REAL :: ustr  (nx,ny,nz)   ! u*rhostr
  REAL :: vstr  (nx,ny,nz)   ! v*rhostr
  REAL :: wstr  (nx,ny,nz)   ! wcont*rhostr

  REAL :: mapfct2(nx,ny)     ! map factor for scalar **2
  REAL :: s1    (0:nx,0:ny,0:nz)! Time lev. 1 of a scalar field to be advected.
  REAL :: s2    (0:nx,0:ny,0:nz)! Time lev. 2 of a scalar field to be advected.
  REAL :: rhostr(nx,ny,nz)   ! rhobar*j3

  REAL :: sadv  (nx,ny,nz)   ! advection term

!
  REAL :: fluxx (nx,ny,nz)   ! temporary work array
  REAL :: fluxy (nx,ny,nz)   ! temporary work array
  REAL :: fluxz (nx,ny,nz)   ! temporary work array
  REAL :: fout  (nx,ny,nz)   ! temporary work array
  REAL :: bout  (0:nx,0:ny,0:nz)   ! temporary work array

!  REAL :: mp_tem(MAX(nx+1,ny+1)*(nz+1))  ! Temporary message passing array.
  REAL :: mp_tem(nx,ny,nz)         ! Temporary message passing array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: tem, eps
!  REAL :: dtxy,dtxz,dtyz,
  REAL :: rdxyz,dxyz

!  INTEGER :: mptag
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'         ! the boundary flags
  INCLUDE 'globcst.inc'       ! Global model control parameters
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  dtxy=dtbig1*dx*dy
!  dtxz=dtbig1*dz*dx
!  dtyz=dtbig1*dz*dy
  dxyz=dx*dy*dz
  rdxyz=1./(dxyz)
  eps=1.e-20
!
!-----------------------------------------------------------------------
!
!  Compute the second-order advective fluxes d(rhostr*u*s)/dx
!
!-----------------------------------------------------------------------
!
  CALL advflxs(nx,ny,nz,nt,dx,dy,dz,dtbig1,                             &
               ustr,vstr,wstr,s2,fctorderopt,                           &
               fluxx,fluxy,fluxz,                                       &
               bout,fout)
! Attention: bout & fout are used as a work arrays here!

!  compute the total outgoing flux for each grid cell (FOUT)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        fout(i,j,k)=rdxyz*( mapfct2(i,j)*(                              &
                           (ABS(fluxx(i  ,j,k))-fluxx(i  ,j,k)) +       &
                           (ABS(fluxx(i+1,j,k))+fluxx(i+1,j,k)) +       &
                           (ABS(fluxy(i  ,j,k))-fluxy(i  ,j,k)) +       &
                           (ABS(fluxy(i,j+1,k))+fluxy(i,j+1,k)))+       &
                           (ABS(fluxz(i,j  ,k))-fluxz(i,j  ,k)) +       &
                           (ABS(fluxz(i,j,k+1))+fluxz(i,j,k+1)))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        s1(i,j,k) = MAX(0.0, s1(i,j,k)) ! Set negative values to zero

        IF(fout(i,j,k) <= s1(i,j,k)*rhostr(i,j,k)) THEN
          bout(i,j,k) = 1.0
        ELSE
          bout(i,j,k) = (s1(i,j,k)*rhostr(i,j,k))/(fout(i,j,k)+eps)
        END IF

      END DO
    END DO
  END DO

!  IF (mp_opt > 0) THEN
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendextew(bout,nx,ny,nz,ebc,wbc,mptag,mp_tem)
!    CALL mprecvextew(bout,nx,ny,nz,ebc,wbc,mptag,mp_tem)
!    CALL mpsendextns(bout,nx,ny,nz,nbc,sbc,mptag,mp_tem)
!    CALL mprecvextns(bout,nx,ny,nz,nbc,sbc,mptag,mp_tem)
!  END IF
  CALL acct_interrupt(bc_acct)
  CALL extndsbc(bout,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)
  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
!  Application of the flux limiter to fluxes
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        IF(fluxx(i,j,k) > 0.0) THEN
          fluxx(i,j,k)=bout(i-1,j,k)*fluxx(i,j,k)
        ELSE
          fluxx(i,j,k)=bout(i,j,k)  *fluxx(i,j,k)
        END IF
      END DO
    END DO
  END DO
!
  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        IF(fluxy(i,j,k) > 0.0) THEN
          fluxy(i,j,k)=bout(i,j-1,k)*fluxy(i,j,k)
        ELSE
          fluxy(i,j,k)=bout(i,j,k)  *fluxy(i,j,k)
        END IF
      END DO
    END DO
  END DO
  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx-1
        IF(fluxz(i,j,k) > 0.0) THEN
          fluxz(i,j,k)=bout(i,j,k-1)*fluxz(i,j,k)
        ELSE
          fluxz(i,j,k)=bout(i,j,k)  *fluxz(i,j,k)
        END IF
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate advection term (sadv)
!
!-----------------------------------------------------------------------

  tem = 1.0/(dtbig1*dx*dy*dz)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        sadv(i,j,k)=sadv(i,j,k)+                                        &
                    ((fluxx(i+1,j,k)-fluxx(i,j,k) +                     &
                      fluxy(i,j+1,k)-fluxy(i,j,k))*mapfct2(i,j) +       &
                      fluxz(i,j,k+1)-fluxz(i,j,k))*tem
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE advmpdcd
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE FCTLIM                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE fctlim(nx,ny,nz,dxyz,std,sa,rhostr,mapfct2,                  &
           fluxx,fluxy,fluxz,cx,cy,cz,rminus,rplus)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform 3-D flux limiting (Zalesak, 1979) on fluxx, fluxy, fluxz
!  which are fluxes in x, y, and z direction respectivly.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/1/1988
!
!  (Originally written for the sigma-coordinate model of
!  Xue and Thorpe, 1991).
!
!  MODIFICATION HISTORY
!
!  11/20/92  (Yifeng Tang)
!  Adapted for ARPS, taking into account of various boundary condition
!  options.
!
!  10/12/1996 (M. Xue)
!  Cleaned up and portions rewritten for ARPS.
!
!  11/20/1997 (Fanyou Kong - CMRP)
!  Added map factor
!
!  7/16/1998 (M. Xue and R. Richardson)
!  Significant changes made. The extrema in the advected
!  (instead of the previous rhostr multipled) field are now used
!  in determining the flux correction coefficients.
!  A one-dimensional prelimiter is included as an option,
!  although it's not recommended for common use.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dxyz     dx*dy*dz
!
!    std      the scalar at n+1 computed from donnor cell fluxes
!    sa       the scalar at time n
!             std and sa are unchanged on return.
!    rhostr   rhobar * j3
!    mapfct2  Map factor for scalars **2
!    fluxx    the antidiffusive flux in x direction
!    fluxy    the antidiffusive flux in y direction
!    fluxz    the antidiffusive flux in z direction
!    cx       the donnor cell flux in x direction
!    cy       the donnor cell flux in y direction
!    cz       the donnor cell flux in z direction
!
!  OUTPUT:
!
!    fluxx    the limited antidiffusive flux in x direction
!    fluxy    the limited antidiffusive flux at y
!    fluxz    the limited antidiffusive flux at z
!
!  WORK ARRAYS:
!
!    rminus    the work array
!    rplus     the work array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: std(0:nx,0:ny,0:nz)
  REAL :: sa (0:nx,0:ny,0:nz)

  REAL :: rhostr(nx,ny,nz)
  REAL :: mapfct2(nx,ny)
  REAL :: fluxx (nx,ny,nz)
  REAL :: fluxy (nx,ny,nz)
  REAL :: fluxz (nx,ny,nz)
  REAL :: cx    (nx,ny,nz)
  REAL :: cy    (nx,ny,nz)
  REAL :: cz    (nx,ny,nz)
  REAL :: rminus(nx,ny,nz)
  REAL :: rplus (nx,ny,nz)
  REAL :: dxyz
  REAL :: qplus,qminus,pplus,pminus

  INTEGER :: fct1dlim  ! Switch for 1D prelimiter
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
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

  INCLUDE 'bndry.inc'           ! the boundary flags
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  fct1dlim = 0

  IF( fct1dlim == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Pre-limiting in x-direciton.
!
!-----------------------------------------------------------------------
!

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          pplus=(MAX(0.,fluxx(i,j,k))-MIN(0.,fluxx(i+1,j,k)))           &
                 *mapfct2(i,j)
          qplus=( MAX(std(i-1,j,k),std(i+1,j,k),std(i,j,k),             &
                      sa (i-1,j,k),sa (i+1,j,k),sa (i,j,k))             &
                      -std(i,j,k) )*dxyz * rhostr(i,j,k)
          IF( pplus < 1.0E-20) THEN
            rplus(i,j,k)=0.0
          ELSE
            rplus(i,j,k)=MIN(1.,qplus/pplus )
          END IF

        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          pminus=(MAX(0.,fluxx(i+1,j,k))-MIN(0.,fluxx(i,j,k)))          &
                  *mapfct2(i,j)
          qminus=-( MIN(std(i-1,j,k),std(i+1,j,k),std(i,j,k),           &
                        sa (i-1,j,k),sa (i+1,j,k),sa (i,j,k))           &
                       -std(i,j,k) )*dxyz * rhostr(i,j,k)

          IF( pminus < 1.e-20 ) THEN
            rminus(i,j,k)=0.0
          ELSE
            rminus(i,j,k)= MIN(1.,qminus/pminus)
          END IF

        END DO
      END DO
    END DO
!
!  Determine the correction coefficient C
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=2,nx-1
          IF( fluxx(i,j,k) >= 0.0) THEN
            cx(i,j,k)=MIN(rplus (i,j,k),rminus(i-1,j,k))
          ELSE
            cx(i,j,k)=MIN(rminus(i,j,k),rplus (i-1,j,k))
          END IF
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Pre-limiting in y-direciton.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          pplus=(MAX(0.,fluxy(i,j,k))-MIN(0.,fluxy(i,j+1,k)))           &
                *mapfct2(i,j)
          qplus=( MAX(std(i,j,k),std(i,j-1,k),std(i,j+1,k),             &
                      sa (i,j,k),sa (i,j-1,k),sa (i,j+1,k))             &
                     -std(i,j,k) )*dxyz * rhostr(i,j,k)
          IF( pplus < 1.0E-20) THEN
            rplus(i,j,k)=0.0
          ELSE
            rplus(i,j,k)=MIN(1.,qplus/pplus )
          END IF
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          pminus=(MAX(0.,fluxy(i,j+1,k))-MIN(0.,fluxy(i,j,k)))          &
                  *mapfct2(i,j)
          qminus=-( MIN(std(i,j,k),std(i,j-1,k),std(i,j+1,k),           &
                      sa (i,j,k),sa (i,j-1,k),sa (i,j+1,k))             &
                      -std(i,j,k) )*dxyz * rhostr(i,j,k)
          IF( pminus < 1.e-20 ) THEN
            rminus(i,j,k)=0.0
          ELSE
            rminus(i,j,k)= MIN(1.,qminus/pminus)
          END IF
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=2,ny-1
        DO i=1,nx-1
          IF( fluxy(i,j,k) >= 0.0) THEN
            cy(i,j,k)=MIN(rplus (i,j,k),rminus(i,j-1,k))
          ELSE
            cy(i,j,k)=MIN(rminus(i,j,k),rplus (i,j-1,k))
          END IF
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Pre-limiting in z-direciton.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          pplus=(MAX(0.,fluxz(i,j,k))-MIN(0.,fluxz(i,j,k+1)))
          qplus=( MAX(std(i,j,k),std(i,j,k-1),std(i,j,k+1),             &
                      sa (i,j,k),sa (i,j,k-1),sa (i,j,k+1))             &
                     -std(i,j,k) )*dxyz * rhostr(i,j,k)
          IF( pplus < 1.0E-20) THEN
            rplus(i,j,k)=0.0
          ELSE
            rplus(i,j,k)=MIN(1.,qplus/pplus )
          END IF
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          pminus=(MAX(0.,fluxz(i,j,k+1))-MIN(0.,fluxz(i,j,k)))
          qminus=-( MIN(std(i,j,k),std(i,j,k-1),std(i,j,k+1),           &
                      sa (i,j,k),sa (i,j,k-1),sa (i,j,k+1))             &
                      -std(i,j,k) )*dxyz * rhostr(i,j,k)

          IF( pminus < 1.e-20 ) THEN
            rminus(i,j,k)=0.0
          ELSE
            rminus(i,j,k)= MIN(1.,qminus/pminus)
          END IF

        END DO
      END DO
    END DO

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF( fluxz(i,j,k) >= 0.0) THEN
            cz(i,j,k)=MIN(rplus (i,j,k),rminus(i,j,k-1))
          ELSE
            cz(i,j,k)=MIN(rminus(i,j,k),rplus (i,j,k-1))
          END IF
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Apply pre-limiting factor to the antidiffusive fluxes.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=2,nx-1
          fluxx(i,j,k)=fluxx(i,j,k)*cx(i,j,k)
        END DO
      END DO
    END DO
!
    DO k=1,nz-1
      DO j=2,ny-1
        DO i=1,nx-1
          fluxy(i,j,k)=fluxy(i,j,k)*cy(i,j,k)
        END DO
      END DO
    END DO
!
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          fluxz(i,j,k)=fluxz(i,j,k)*cz(i,j,k)
        END DO
      END DO
    END DO

  END IF  ! End of 1D prelimiting step

!-----------------------------------------------------------------------
!
!  Calculate total flux directing towards point (i,j,k) pplus
!  then inward flux amplification factor r(i,j,k)+
!  r+=min(1.,q+/p+).
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        pplus=(MAX(0.,fluxx(i,j,k))-MIN(0.,fluxx(i+1,j,k))              &
              +MAX(0.,fluxy(i,j,k))-MIN(0.,fluxy(i,j+1,k)))             &
              *mapfct2(i,j)                                             &
              +MAX(0.,fluxz(i,j,k))-MIN(0.,fluxz(i,j,k+1))
        qplus=(MAX(std(i-1,j,k),std(i+1,j,k),std(i,j,k),                &
                   std(i,j-1,k),std(i,j+1,k),std(i,j,k-1),std(i,j,k+1), &
                   sa (i-1,j,k),sa (i+1,j,k),sa (i,j,k),                &
                   sa (i,j-1,k),sa (i,j+1,k),sa (i,j,k-1),sa (i,j,k+1)) &
              -std(i,j,k) )*dxyz * rhostr(i,j,k)
        IF( pplus < 1.0E-20) THEN
          rplus(i,j,k)=0.0
        ELSE
          rplus(i,j,k)=MIN(1.,qplus/pplus )
        END IF

      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate total flux directing away from point (i,j,k) pminus
!  then outward flux amplification factor r(i,j,k)-
!  r- = min(1.,q-/p-).
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        pminus=(MAX(0.,fluxx(i+1,j,k))-MIN(0.,fluxx(i,j,k))             &
               +MAX(0.,fluxy(i,j+1,k))-MIN(0.,fluxy(i,j,k)))            &
               *mapfct2(i,j)                                            &
               +MAX(0.,fluxz(i,j,k+1))-MIN(0.,fluxz(i,j,k))
        qminus=-(MIN(std(i-1,j,k),std(i+1,j,k),std(i,j,k),              &
                    std(i,j-1,k),std(i,j+1,k),std(i,j,k-1),std(i,j,k+1), &
                    sa (i-1,j,k),sa (i+1,j,k),sa (i,j,k),               &
                    sa (i,j-1,k),sa (i,j+1,k),sa (i,j,k-1),sa (i,j,k+1)) &
                    -std(i,j,k) )*dxyz * rhostr(i,j,k)

        IF( pminus < 1.e-20 ) THEN
          rminus(i,j,k)=0.0
        ELSE
          rminus(i,j,k)= MIN(1.,qminus/pminus)
        END IF
      END DO
    END DO
  END DO
!
!  Determine the correction coefficient C
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=2,nx-1
        IF( fluxx(i,j,k) >= 0.0) THEN
          cx(i,j,k)=MIN(rplus (i,j,k),rminus(i-1,j,k))
        ELSE
          cx(i,j,k)=MIN(rminus(i,j,k),rplus (i-1,j,k))
        END IF
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=2,ny-1
      DO i=1,nx-1
        IF( fluxy(i,j,k) >= 0.0) THEN
          cy(i,j,k)=MIN(rplus (i,j,k),rminus(i,j-1,k))
        ELSE
          cy(i,j,k)=MIN(rminus(i,j,k),rplus (i,j-1,k))
        END IF
      END DO
    END DO
  END DO

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        IF( fluxz(i,j,k) >= 0.0) THEN
          cz(i,j,k)=MIN(rplus (i,j,k),rminus(i,j,k-1))
        ELSE
          cz(i,j,k)=MIN(rminus(i,j,k),rplus (i,j,k-1))
        END IF
      END DO
    END DO
  END DO
!
!  Apply final correction factors to the antidiffusive fluxes
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=2,nx-1
        fluxx(i,j,k)=fluxx(i,j,k)*cx(i,j,k)
      END DO
    END DO
  END DO
!
  DO k=1,nz-1
    DO j=2,ny-1
      DO i=1,nx-1
        fluxy(i,j,k)=fluxy(i,j,k)*cy(i,j,k)
      END DO
    END DO
  END DO
!
  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        fluxz(i,j,k)=fluxz(i,j,k)*cz(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE fctlim
