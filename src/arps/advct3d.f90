!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ADVUVW                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE advuvw(nx,ny,nz,u,v,w,wcont,rhostr,ubar,vbar, mapfct,        &
           uforce,vforce,wforce,uadv,vadv,wadv,                         &
           ustr,vstr,wstr,tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinates the calculation of the advection terms uadv, vadv and
!  wadv of the u, v and w equations. These terms are written in
!  equivalent advection form.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (M. Xue)
!  Added full documentation.
!
!  5/29/92 (K. Brewster)
!  Further facelift.
!
!  4/9/93 (M. Xue & K. Brewster)
!  Some index bounds for operators and loop bounds were corrected.
!  Redundant calculations were done before at the boundaries that
!  caused out-of-bound calculations, which should not have
!  affected the results in normal situations though.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  2/9/94 (D. Jahn)
!  Add tem3 work array in conjunction w/ changes to advu,advv,
!  advw, and advcts.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/11/98 (D. Weber)
!  Removed most single operation loops and merged the advection
!  contributions into the forcing arrays prior to the return.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mapfct   Map factors at scalar, u and v points
!
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates at all time levels (m/s).
!    wcont    Contravariant vertical velocity (m/s)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    ubar     Base state x component of velocity (m/s)
!    vbar     Base state y component of velocity (m/s)
!
!  OUTPUT:
!
!    uadv     u-equation advection terms  (kg/m**3)*(m/s)/s
!    vadv     v-equation advection terms  (kg/m**3)*(m/s)/s
!    wadv     w-equation advection terms  (kg/m**3)*(m/s)/s
!
!  WORK ARRAYS:
!
!    ustr     u * rhostr
!    vstr     v * rhostr
!    wstr     w * rhostr
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

  INCLUDE 'timelvls.inc'

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)

  REAL :: uadv  (nx,ny,nz)     ! u-eqn advection terms
                               ! (kg/m**3)*(m/s)/s
  REAL :: vadv  (nx,ny,nz)     ! v-eqn advection terms
                               ! (kg/m**3)*(m/s)/s
  REAL :: wadv  (nx,ny,nz)     ! w-eqn advection terms
                               ! (kg/m**3)*(m/s)/s

  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)
                               ! uforce= -uadv + umix + ucorio

  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)
                               ! vforce= -vadv + vmix + vcorio

  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
                               ! wforce= -wadv + wmix + wbuoy

  REAL :: ustr  (nx,ny,nz)     ! u * rhostr
  REAL :: vstr  (nx,ny,nz)     ! v * rhostr
  REAL :: wstr  (nx,ny,nz)     ! w * rhostr

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
!  REAL :: tema
  INTEGER :: tlevel
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tlevel=tpresent
!
!-----------------------------------------------------------------------
!
!  Calculate ustr=rhostr*u, vstr=rhostr*v, wstr=rhostr*w
!
!-----------------------------------------------------------------------
!
  CALL uvwrho(nx,ny,nz,u(1,1,1,tlevel),v(1,1,1,tlevel),                 &
              wcont,rhostr,                                             &
              ustr,vstr,wstr)
!
!-----------------------------------------------------------------------
!
!  Calculate uadv:
!
!-----------------------------------------------------------------------
!
  CALL advu(nx,ny,nz, u,v,w,                                            &
            ustr,vstr,wstr, ubar, mapfct(1,1,2), uadv, uforce,          &
            tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  Calculate vadv:
!
!-----------------------------------------------------------------------
!
  CALL advv(nx,ny,nz, u,v,w,                                            &
            ustr,vstr,wstr, vbar, mapfct(1,1,3), vadv, vforce,          &
            tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  Calculate wadv:
!
!-----------------------------------------------------------------------
!
  CALL advw(nx,ny,nz, u,v,w,                                            &
            ustr,vstr,wstr, mapfct(1,1,1), wadv, wforce,                &
            tem1,tem2,tem3)

  RETURN
END SUBROUTINE advuvw
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ADVP                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advp(nx,ny,nz,pprt,u,v,w,wcont,rhostr, mapfct,               &
           j3,aj3x,aj3y,aj3z,                                           &
           padv,                                                        &
           uj3,vj3,wj3,tem1,tem2,tem3,tem4,mp_tem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the advection of perturbation pressure.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (M. Xue)
!  Added full documentation.
!
!  5/29/92 (K. Brewster)
!  Further facelift.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  5/25/93 (M. Xue)
!  Fixed the bug in the base state vertical pressure term.
!  A j3 was missing before.
!
!  5/26/93 (M. Xue)
!  Corrected pbar advection term.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  9/7/94 (M. Xue)
!  This subroutine now calculates the perturbation pressure advection
!  only, by calling subroutine ADVCTS.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  9/18/98 (D. Weber)
!  Added aj3x,y,z arrays.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    pprt     Perturbation pressure at all time levels (Pascal)
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates at al time levels (m/s).
!    wcont    Contravariant vertical velocity (m/s)
!
!    rhostr   Base state air density (kg/m**3)
!
!    mapfct   Map factors at scalar points
!
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    padv     Pressure equation advection terms (kg/(m*s**3))
!
!  WORK ARRAYS:
!
!    uj3      Temporary work array.
!    vj3      Temporary work array.
!    wj3      Temporary work array.
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
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

  INCLUDE 'timelvls.inc'

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure from that
                               ! of base state atmosphere (Pascal)

  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3)

  REAL :: mapfct(nx,ny)        ! Map factors at scalar points

  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.

  REAL :: padv  (nx,ny,nz)     ! Advection term of the pressure
                               ! equation (kg/(m*s**3)

  REAL :: uj3 (nx,ny,nz)       ! Temporary work array to carry u*j3
  REAL :: vj3 (nx,ny,nz)       ! Temporary work array to carry v*j3
  REAL :: wj3 (nx,ny,nz)       ! Temporary work array to carry wcont*j3
  REAL :: tem1(nx,ny,nz)       ! Temporary work array
  REAL :: tem2(nx,ny,nz)       ! Temporary work array
  REAL :: tem3(nx,ny,nz)       ! Temporary work array
  REAL :: tem4(nx,ny,nz)       ! Temporary work array

  REAL :: mp_tem(nx,ny,nz)     ! Temporary message passing array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: tlevel
!  REAL :: tema
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
  tlevel=tpresent
!
!-----------------------------------------------------------------------
!
!  Advection of the perturbation pressure
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Calculate the advection of perturbation pressure.
!  (j3*u*dp/dx+j3*v*dp/dy+j3*wcont*dp/dz) =
!    avgx( avgx(j3) * u * difx(p) )
!  + avgy( avgy(j3) * u * dify(p) )
!  + avgz( avgz(j3) * wcont * difz(p) )
!
!  Subroutine ADVCTS is called to evaluate this term.
!
!-----------------------------------------------------------------------


  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        uj3(i,j,k)=u(i,j,k,tlevel)*aj3x(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        vj3(i,j,k)=v(i,j,k,tlevel)*aj3y(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx-1
        wj3(i,j,k)=wcont(i,j,k)*aj3z(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Perturbation potential temperature advection
!
!-----------------------------------------------------------------------
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem4(i,j,k) = 0.0
      END DO
    END DO
  END DO

  CALL advcts(nx,ny,nz, pprt,                                           &
              u(1,1,1,tlevel),v(1,1,1,tlevel),                          &
              uj3,vj3,wj3,tem4, mapfct,                                 &
              padv, tem1,tem2,tem3,mp_tem)
!
!-----------------------------------------------------------------------
!
!  Note that the base state pressure advection term is evaluated
!  inside the small time step, therefore padv does not include
!  this part on exit of this subroutine.
!
!-----------------------------------------------------------------------
!
  RETURN
END SUBROUTINE advp

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADVQ                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advq(nx,ny,nz,q,u,v,wcont,ustr,vstr,wstr,                    &
           rhostr,qbar, mapfct,                                         &
           qadv,                                                        &
           tem1,tem2,tem3,mp_tem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate water/ice advection.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  05/20/92.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (M. Xue)
!  Added full documentation.
!
!  5/29/92 (K. Brewster)
!  Further facelift and order of argument list changed.
!
!  2/9/94 (D. Jahn)
!  Add tem3 to the call of advcts.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  9/28/98 (Dan Weber)
!  Moved u,v,wstr computations outside of this routine.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mapfct   Map factors at scalar points
!
!    q        Water or ice mixing ratio (kg/kg)
!
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    rhostr   Base state air density (kg/m**3)
!    qbar     Base state of water oy ice mixing ratio (kg/kg)
!
!  OUTPUT:
!
!    qadv     Water or ice mixing ratio advection term in the
!             conservation equation (kg/m**3)*(kg/kg)/s
!
!  WORK ARRAYS:
!
!    ustr     u * rhostr
!    vstr     v * rhostr
!    wstr     w * rhostr
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
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

  INCLUDE 'timelvls.inc'

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: mapfct(nx,ny)        ! Map factors at scalar points

  REAL :: q     (nx,ny,nz,nt)  ! One of the water/ice variables (kg/kg)
  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)

  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3)

  REAL :: qbar  (nx,ny,nz)     ! Base state of the water/ice variable
  REAL :: qadv  (nx,ny,nz)     ! Advection of a water/ice variable

  REAL :: ustr  (nx,ny,nz)     ! u * rhostr
  REAL :: vstr  (nx,ny,nz)     ! v * rhostr
  REAL :: wstr  (nx,ny,nz)     ! w * rhostr

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array

  REAL :: mp_tem(nx,ny,nz)     ! Temporary message passing array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: tlevel
!  INTEGER :: i, j, k
!  REAL :: tema
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global model control parameters
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tlevel=tpresent
!
!-----------------------------------------------------------------------
!
!  Calculate ustr=rhostr*u, vstr=rhostr*v, wstr=rhostr*w
!
!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!
!  Advection of the water/ice mixing ratio
!
!-----------------------------------------------------------------------

  CALL advcts(nx,ny,nz, q,                                              &
              u(1,1,1,tlevel),v(1,1,1,tlevel),                          &
              ustr,vstr,wstr, qbar, mapfct,                             &
              qadv, tem1,tem2,tem3,mp_tem)

  RETURN
END SUBROUTINE advq
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADVU                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advu(nx,ny,nz,u,v,w,ustr,vstr,wstr, ubar, mapfct,            &
           uadv, uforce,                                                &
           tem1,tem2,tem3)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the advection terms of the u-equation. These terms are
!  written in equivalent advection form.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (M. Xue)
!  Added full documentation.
!
!  5/29/92 (K. Brewster)
!  Further facelift.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  2/9/94 (D. Jahn)
!  Add 4th-order momentum advection.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  11/5/1995 (M. Xue)
!  Added vertical fourth order advection for u.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  9/1/98 (D. Weber and T. Chung)
!  Removed operators and merged loops to improve code efficiency.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mapfct   Map factors at u points
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at all time levels (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates at all time levels (m/s).
!
!    ustr     u * rhostr
!    vstr     v * rhostr
!    wstr     w * rhostr
!
!    ubar     Base state x component of velocity (m/s)
!
!  OUTPUT:
!
!    uforce   Acoustically inactive forcing terms in u-momentum
!             equation (kg/(m*s)**2).
!
!  WORK ARRAYS:
!
!    uadv     u-equation advection terms (kg/m**3)*(m/s)/s
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

  INCLUDE 'timelvls.inc'

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: mapfct(nx,ny)        ! Map factors at u points

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity at a given time level
                               ! (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity at a given time level
                               ! (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity at a given time level
                               ! (m/s)
  REAL :: ustr  (nx,ny,nz)     ! u*rhostr
  REAL :: vstr  (nx,ny,nz)     ! v*rhostr
  REAL :: wstr  (nx,ny,nz)     ! w*rhostr
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)

  REAL :: uadv  (nx,ny,nz)     ! Advection term of the u-eq.
                               ! (kg/m**3)*(m/s)/s

  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: ibgn,iend,jbgn,jend
  REAL :: vsb, vnb
  REAL :: fourthirds,foursixth,onesix
  REAL :: tema

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'       ! Global model control parameters
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'bndry.inc'         ! Parameters for boundary conditions.
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  fourthirds = 4.0/3.0
  foursixth = 4.0/6.0
  onesix = 1.0/6.0

!-----------------------------------------------------------------------
!
!  Calculate rho*u* du/dx = avgx( avgx( ustr )*difx( u ) ):
!
!-----------------------------------------------------------------------

  tema = 0.25*dxinv
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=2,nx-1
        uadv(i,j,k)=tema*((ustr(i+1,j,k)+ustr(i,j,k))*                  &
                          (u(i+1,j,k,2)-u(i,j,k,2))                     &
                         +(ustr(i,j,k)+ustr(i-1,j,k))*                  &
                          (u(i,j,k,2)-u(i-1,j,k,2)))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!   tem2 contains avgx((avgx(u*rho))*du/dx)  (second order term)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  If momentum advection 4th-order, calculate adjustment term:
!  avg2x(avg2x(ustr)*dif2x(u))
!
!-----------------------------------------------------------------------

  IF (madvopt == 2 .OR. madvopt == 3) THEN

    tema = 0.25*dxinv
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-1
          tem2(i,j,k)=tema*(ustr(i+1,j,k)+ustr(i-1,j,k))*               &
                            (u(i+1,j,k,2)-u(i-1,j,k,2))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!   At this point tem2 contains avg2x(ustr)*dif2x(u)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!   In the case of periodic boundary condition, need to calculate
!   the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem2,nx,ny,nz,ebc,wbc,1,tem1)
    END IF

    CALL acct_interrupt(bc_acct)

    ibgn = 3
    iend = nx-2
    jbgn = 1
    jend = ny-1

    IF (wbc == 0) ibgn=2

    IF (wbc == 2) THEN
      ibgn=2
      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO j=jbgn,jend
          tem2(1,j,k) = tem2(nx-2,j,k)
        END DO
        END DO
      END IF
    END IF

    IF (ebc == 0) iend = nx-1

    IF (ebc == 2) THEN
      iend = nx-1
      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO j=jbgn,jend
          tem2(nx,j,k) = tem2(3,j,k)
        END DO
        END DO
      END IF
    END IF


!-----------------------------------------------------------------------
!
!  In the case of rigid boundary condition, need to calculate
!  the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (wbc == 1) THEN
      ibgn=2
      DO k=2,nz-2
        DO j=jbgn,jend
          tem2(1,j,k) = -tem2(3,j,k)
        END DO
      END DO
    END IF

    IF (ebc == 1) THEN
      iend = nx-1
      DO k=2,nz-2
        DO j=jbgn,jend
          tem2(nx,j,k) = -tem2(nx-2,j,k)
        END DO
      END DO
    END IF

    CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Concatenate 2nd-order and adjusting term for 4th-order
!
!-----------------------------------------------------------------------

    DO k=2,nz-2
      DO j=jbgn,jend
        DO i=ibgn,iend
          uadv(i,j,k) = fourthirds*uadv(i,j,k) -                        &
                        onesix*(tem2(i+1,j,k)+tem2(i-1,j,k))
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Calculate rho*v* du/dy = avgy( avgx( vstr )* dify( u ) )
!
!-----------------------------------------------------------------------

  tema = 0.25*dyinv
  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-1
        tem2(i,j,k)=tema*((vstr(i,j,k)+vstr(i-1,j,k))*                  &
                          (u(i,j,k,2)-u(i,j-1,k,2))                     &
                         +(vstr(i,j+1,k)+vstr(i-1,j+1,k))*              &
                          (u(i,j+1,k,2)-u(i,j,k,2)))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  tem1 contains avgy( avgx( vstr )* dify( u ) )  second order term
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!    For 4th-order momentum advection, calculate the adjustment term.
!        avg2y( avgx( avgy (vstr) )* dif2y( u ))
!
!-----------------------------------------------------------------------

  IF (madvopt == 2 .OR. madvopt == 3) THEN

    DO k=2,nz-2     ! 2 bugs found here..multiply sb plus...
      DO j=2,ny-2
        DO i=2,nx-1
          tem1(i,j,k)=0.25*((vstr(i,j+1,k)+vstr(i,j,k))                 &
                           +(vstr(i-1,j+1,k)+vstr(i-1,j,k)))
        END DO
      END DO
    END DO

    tema = 0.5*dyinv
    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-1
          tem3(i,j,k)=tema*(u(i,j+1,k,2)-u(i,j-1,k,2))*tem1(i,j,k)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  At this point tem3 contains avgx( avgy (vstr) )* dif2y( u ))
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!    In the case of periodic boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dns(tem3,nx,ny,nz,nbc,sbc,0,tem1)
    END IF

    CALL acct_interrupt(bc_acct)

    ibgn = 2
    iend = nx-1
    jbgn = 3
    jend = ny-3

    IF (nbc == 0) jend = ny-2

    IF (nbc == 2) THEN
      jend = ny-2
      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,ny-1,k) = tem3(i,2,k)
        END DO
        END DO
      END IF
    END IF

    IF (sbc == 0) jbgn = 2

    IF (sbc == 2) THEN
      jbgn = 2
      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,1,k) = tem3(i,ny-2,k)
        END DO
        END DO
      END IF
    END IF

!-----------------------------------------------------------------------
!
!    In the case of rigid boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (nbc == 1) THEN
      jend = ny-2
      DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,ny-1,k) = tem3(i,ny-2,k)
        END DO
      END DO
    END IF

    IF (sbc == 1) THEN
      jbgn = 2
      DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,1,k) = tem3(i,2,k)
        END DO
      END DO
    END IF

    CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  concatenate 2nd-order and adjusting term for 4th-order
!
!-----------------------------------------------------------------------

    DO k=2,nz-2
      DO j=jbgn,jend
        DO i=ibgn,iend
          tem2(i,j,k) = fourthirds*tem2(i,j,k) -                        &
                       onesix*(tem3(i,j+1,k)+tem3(i,j-1,k))
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Calculate rho*v* du/dy on the north and south boundaries using
!  one sided advection.
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO i=2,nx-1
      vsb=(vstr(i-1,1,k)+vstr(i-1,2,k)                                  &
          +vstr(i,1,k)+vstr(i,2,k))*0.25
      vnb=(vstr(i-1,ny-1,k)+vstr(i-1,ny,k)                              &
          +vstr(i,ny-1,k)+vstr(i,ny,k))*0.25
      IF (vsb < 0.0) THEN
        tem2(i,1,k)=vsb*(u(i,2,k,tpresent)-u(i,1,k,tpresent))*dyinv
      ELSE
        tem2(i,1,k)= rlxlbc* vsb*(u(i,1,k,tpast)-ubar(i,1,k))*dyinv
      END IF

      IF (vnb > 0.0) THEN
        tem2(i,ny-1,k)=vnb*                                             &
               (u(i,ny-1,k,tpresent)-u(i,ny-2,k,tpresent))*dyinv
      ELSE
        tem2(i,ny-1,k)=-rlxlbc*vnb*                                     &
               (u(i,ny-1,k,tpast)-ubar(i,ny-1,k))*dyinv
      END IF
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Add the u-equation horizontal advection terms, rho*u*du/dx (uadv) and
!  rho*v*du/dy (tem1).
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=2,nx-1
        uforce(i,j,k)=uforce(i,j,k) -                                   &
                       (uadv(i,j,k)+tem2(i,j,k)) * mapfct(i,j)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate rho*w* du/dz = avgz( avgx( wstr ) * difz ( u ) )
!
!-----------------------------------------------------------------------

  tema=0.5*dzinv
  DO k=2,nz-1
    DO j=1,ny-1
      DO i=2,nx-1
        tem2(i,j,k)=tema*(wstr(i,j,k)+wstr(i-1,j,k))*                   &
                          (u(i,j,k,2)-u(i,j,k-1,2))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  tem2 contains avgz( avgx( wstr ) * difz ( u ) )
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  For 4th-order momentum advection, calculate the fourth order
!  contribution  avg2z( avgz( avgx (wstr) )* dif2z( u ))
!
!-----------------------------------------------------------------------

  IF (madvopt == 3) THEN

    tema = 0.125*dzinv
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-1
          tem3(i,j,k)=tema*((wstr(i,j,k+1)+wstr(i-1,j,k+1))+            &
                           (wstr(i,j,k)+wstr(i-1,j,k))) *               &
                           (u(i,j,k+1,2)-u(i,j,k-1,2))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Add the fourth order adjustment term to the 2nd order term
!  and add the u-equation vertical advection term to the
!  horizontal advection terms.
!
!-----------------------------------------------------------------------

    DO k=3,nz-3
      DO j=1,ny-1
        DO i=2,nx-1
          uforce(i,j,k) = uforce(i,j,k) -                               &
                       foursixth*(tem2(i,j,k)+tem2(i,j,k+1)) +          &
                         onesix*(tem3(i,j,k+1)+tem3(i,j,k-1))
        END DO
      END DO
    END DO

    DO k=2,nz-2,nz-4
      DO j=1,ny-1
        DO i=2,nx-1
          uforce(i,j,k)=uforce(i,j,k)-(tem2(i,j,k+1)+tem2(i,j,k))*0.5
        END DO
      END DO
    END DO

  ELSE   !  perform second order advection....

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-1
          uforce(i,j,k)=uforce(i,j,k)-(tem2(i,j,k+1)+tem2(i,j,k))*0.5
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE advu


!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ADVV                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advv(nx,ny,nz,u,v,w,ustr,vstr,wstr, vbar, mapfct,            &
           vadv, vforce,                                                &
           tem1,tem2,tem3)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the advection terms of the v-equation. These terms are
!  written in equivalent advection form.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (M. Xue)
!  Added full documentation.
!
!  5/29/92 (K. Brewster)
!  Further facelift.
!
!  4/9/93 (M. Xue & K. Brewster)
!  Some index bounds for operators and loop bounds were corrected.
!  Redundant calculations were done before at the boundaries that
!  caused out-of-bound calculations, which should not have
!  affected the results in normal situations though.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  2/9/94 (D. Jahn)
!  Add 4th-order momentum advection.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  11/5/1995 (M. Xue)
!  Added vertical fourth order advection for v.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  9/1/98 (D. Weber and T. Chung)
!  Removed operators and merged loops to improve code efficiency.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mapfct   Map factors at v points
!
!    v        y component of velocity at a given time level (m/s)
!
!    ustr     u * rhostr
!    vstr     v * rhostr
!    wstr     w * rhostr
!
!    vbar     Base state y component of velocity  (m/s)
!
!  OUTPUT:
!
!    vforce   Acoustically inactive forcing terms in v-momentum
!             equation (kg/(m*s)**2).
!
!  WORK ARRAYS:
!
!    vadv     v equation advection terms (kg/m**3)*(m/s)/s
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INCLUDE 'timelvls.inc'

  REAL :: mapfct(nx,ny)        ! Map factors at v points

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity at a given time level
                               ! (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity at a given time level
                               ! (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity at a given time level
                               ! (m/s)
  REAL :: ustr  (nx,ny,nz)     ! u*rhostr
  REAL :: vstr  (nx,ny,nz)     ! v*rhostr
  REAL :: wstr  (nx,ny,nz)     ! w*rhostr

  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)

  REAL :: vadv  (nx,ny,nz)     ! Advection term of the v-eq.
                               ! (kg/m**3)*(m/s)/s

  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: ibgn,iend,jbgn,jend
  REAL :: ueb, uwb
  REAL :: fourthirds,foursixth,onesix
  REAL :: tema

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'    ! Global model control parameters
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'bndry.inc'      ! Parameters for boundary conditions.
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  fourthirds = (4.0/3.0)
  foursixth = (4.0/6.0)
  onesix = (1.0/6.0)

!-----------------------------------------------------------------------
!
!  Calculate rho*u* dv/dx  = avgx( avgy( ustr ) * difx( v ) )
!
!-----------------------------------------------------------------------

  tema = 0.25*dxinv
  DO k=2,nz-1
    DO j=2,ny-1
      DO i=2,nx-2
        vadv(i,j,k)=tema*((ustr(i+1,j,k)+ustr(i+1,j-1,k))*              &
                          (v(i+1,j,k,2)-v(i,j,k,2))                     &
                         +(ustr(i,j,k)+ustr(i,j-1,k))*                  &
                          (v(i,j,k,2)-v(i-1,j,k,2)))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!    For 4th-order momentum advection, calculate the adjustment term
!           avg2x( avgx(avgy(ustr))*dif2x(v) )
!
!-----------------------------------------------------------------------

  IF (madvopt == 2 .OR. madvopt == 3) THEN

    tema = 0.125 * dxinv
    DO k=2,nz-2
      DO j=2,ny-1
        DO i=2,nx-2
          tem2(i,j,k)=tema*((ustr(i,j,k)+ustr(i,j-1,k)) +               &
                            (ustr(i+1,j,k)+ustr(i+1,j-1,k))) *          &
                            (v(i+1,j,k,2)-v(i-1,j,k,2))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!    In the case of periodic boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem2,nx,ny,nz,ebc,wbc,0,tem1)
    END IF

    CALL acct_interrupt(bc_acct)

    ibgn = 3
    iend = nx-3
    jbgn = 2
    jend = ny-1

    IF (wbc == 0) ibgn = 2

    IF (wbc == 2) THEN
      ibgn = 2

      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO j=jbgn,jend
          tem2(1,j,k) = tem2(nx-2,j,k)
        END DO
        END DO
      END IF
    END IF

    IF (ebc == 0) iend = nx-2

    IF (ebc == 2) THEN
      iend = nx-2

      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO j=jbgn,jend
          tem2(nx-1,j,k) = tem2(2,j,k)
        END DO
        END DO
      END IF
    END IF

!-----------------------------------------------------------------------
!
!    In the case of rigid boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (wbc == 1) THEN
      ibgn = 2
      DO k=2,nz-2
        DO j=jbgn,jend
          tem2(1,j,k) = tem2(2,j,k)
        END DO
      END DO
    END IF

    IF (ebc == 1) THEN
      iend = nx-2
      DO k=2,nz-2
        DO j=jbgn,jend
          tem2(nx-1,j,k) = tem2(nx-2,j,k)
        END DO
      END DO
    END IF

    CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Concatenate the 2nd and 4th-order terms for urho*dv/dx
!
!-----------------------------------------------------------------------

    DO  k=2,nz-2
      DO  j=jbgn,jend
        DO  i=ibgn,iend
          vadv(i,j,k) = fourthirds*vadv(i,j,k) -                        &
                         onesix*(tem2(i+1,j,k)+tem2(i-1,j,k))
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Calculate rho*u* dv/dx at the east and west boundaries
!  using one-sided advection.
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=2,ny-1
      uwb=(ustr(1,j,k)+ustr(2,j,k)+ustr(1,j-1,k)+ustr(2,j-1,k))*0.25
      ueb=(ustr(nx-1,j,k)+ustr(nx,j,k)+ustr(nx-1,j-1,k)                 &
          +ustr(nx,j-1,k))*0.25
      IF (uwb < 0.0) THEN
        vadv(1,j,k)=uwb*(v(2,j,k,tpresent)-v(1,j,k,tpresent))*dxinv
      ELSE
        vadv(1,j,k)=rlxlbc*uwb*(v(1,j,k,tpast)-vbar(1,j,k))*dxinv
      END IF

      IF (ueb > 0.0) THEN
        vadv(nx-1,j,k)=ueb*                                             &
            (v(nx-1,j,k,tpresent)-v(nx-2,j,k,tpresent))*dxinv
      ELSE
        vadv(nx-1,j,k)=-rlxlbc*ueb*                                     &
            (v(nx-1,j,k,tpast)-vbar(nx-1,j,k))*dxinv
      END IF
    END DO
  END DO


!-----------------------------------------------------------------------
!
!  Calculate rho*v* dv/dy = avgy( avgy( vstr ) * dify( v ) )
!
!-----------------------------------------------------------------------

  tema = 0.25*dyinv
  DO k=2,nz-2
    DO j=2,ny-1
      DO i=1,nx-1
        tem2(i,j,k)=tema*((vstr(i,j+1,k)+vstr(i,j,k))*                  &
                          (v(i,j+1,k,2)-v(i,j,k,2))                     &
                         +(vstr(i,j,k)+vstr(i,j-1,k))*                  &
                          (v(i,j,k,2)-v(i,j-1,k,2)))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  tem2 contains avgy( avgy( vstr ) * dify( v ) ) second order term
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!    For 4th-order momentum advection, calculate the adjustment term:
!            avg2y (avg2y(vstr) * dif2y(v) )
!
!-----------------------------------------------------------------------

  IF (madvopt == 2 .OR. madvopt == 3) THEN

    tema = 0.25*dyinv
    DO k=2,nz-2
      DO j=2,ny-1
        DO i=1,nx-1
          tem3(i,j,k)=tema*(vstr(i,j+1,k)+vstr(i,j-1,k))*               &
                            (v(i,j+1,k,2)-v(i,j-1,k,2))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  At this point tem3 contains avg2y(vstr) * dif2y(v)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!    In the case of periodic boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dns(tem3,nx,ny,nz,nbc,sbc,2,tem1)
    END IF

    CALL acct_interrupt(bc_acct)

    ibgn = 1
    iend = nx-1
    jbgn = 3
    jend = ny-2

    IF (nbc == 0) jend = ny-1

    IF (nbc == 2) THEN
      jend = ny-1
      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,ny,k) = tem3(i,3,k)
        END DO
        END DO
      END IF
    END IF

    IF (sbc == 0) jbgn = 2

    IF (sbc == 2) THEN
      jbgn = 2
      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,1,k) = tem3(i,ny-2,k)
        END DO
        END DO
      END IF
    END IF

!-----------------------------------------------------------------------
!
!    In the case of rigid boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (nbc == 1) THEN
      jend = ny-1
      DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,ny,k) = -tem3(i,ny-2,k)
        END DO
      END DO
    END IF

    IF (sbc == 1) THEN
      jbgn = 2
      DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,1,k) = -tem3(i,3,k)
        END DO
      END DO
    END IF

    CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Concatenate the 2nd- and 4th-order terms:  vstr*dv/dy
!
!-----------------------------------------------------------------------

    DO k=2,nz-2
      DO j=jbgn,jend
        DO i=ibgn,iend
          tem2(i,j,k) = fourthirds*tem2(i,j,k) -                        &
                         onesix*(tem3(i,j+1,k)+tem3(i,j-1,k))
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Sum the x and y contributions to v advection
!  At this point tem1 contains the y contribution to the v advection
!
!-----------------------------------------------------------------------
!
!  Add the v-equation horizontal advection terms, rho*u*dv/dx (vadv) and
!  rho*v*dv/dy (tem1).
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=2,ny-1
      DO i=1,nx-1
        vforce(i,j,k)=vforce(i,j,k) -                                   &
                       (vadv(i,j,k)+tem2(i,j,k)) * mapfct(i,j)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate rho*w* dv/dz = avgz( avgy( wstr ) * difz( v ) )
!
!-----------------------------------------------------------------------

  tema = 0.5*dzinv
  DO k=2,nz-1
    DO j=2,ny-1
      DO i=1,nx-1
        tem2(i,j,k)=tema*(wstr(i,j,k)+wstr(i,j-1,k))*                   &
                          (v(i,j,k,2)-v(i,j,k-1,2))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  tem2 contains avgz( avgy( wstr ) * difz( v ) )
!
!-----------------------------------------------------------------------

  IF (madvopt == 3) THEN

!-----------------------------------------------------------------------
!
!    For 4th-order momentum advection, calculate the adjustment term.
!    avg2z( avgz( avgy (wstr) )* dif2z( v ))
!
!-----------------------------------------------------------------------

    tema = 0.125 * dzinv
    DO k=2,nz-2
      DO j=2,ny-1
        DO i=1,nx-1
          tem3(i,j,k) = tema*((wstr(i,j,k+1)+wstr(i,j-1,k+1))+          &
                              (wstr(i,j,k)+wstr(i,j-1,k))) *            &
                              (v(i,j,k+1,2)-v(i,j,k-1,2))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Add the fourth order contribution to the 2nd order term and
!  add the v-equation vertical advection term rho*w*dv/dz (tem1) to
!  the horizontal advection terms.
!
!-----------------------------------------------------------------------

    DO k=3,nz-3
      DO j=2,ny-1
        DO i=1,nx-1
          vforce(i,j,k) = vforce(i,j,k) -                               &
                      foursixth*(tem2(i,j,k+1)+tem2(i,j,k)) +           &
                         onesix*(tem3(i,j,k+1)+tem3(i,j,k-1))
        END DO
      END DO
    END DO

    DO k=2,nz-2,nz-4
      DO j=2,ny-1
        DO i=1,nx-1
          vforce(i,j,k)=vforce(i,j,k) - (tem2(i,j,k+1)+tem2(i,j,k))*0.5
        END DO
      END DO
    END DO

  ELSE   ! perform second order advection...

    DO k=2,nz-2
      DO j=2,ny-1
        DO i=1,nx-1
          vforce(i,j,k)=vforce(i,j,k) - (tem2(i,j,k+1)+tem2(i,j,k))*0.5
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE advv


!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ADVW                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advw(nx,ny,nz,u,v,w,ustr,vstr,wstr, mapfct,                  &
           wadv, wforce,                                                &
           tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!
!  Calculate the advection terms of the w-equation. These terms are
!  written in equivalent advection form.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (M. Xue)
!  Added full documentation.
!
!  5/29/92 (K. Brewster)
!  Further facelift.
!
!  4/9/93 (M. Xue & K. Brewster)
!  Some index bounds for operators and loop bounds were corrected.
!  Redundant calculations were done before at the boundaries that
!  caused out-of-bound calculations, which should not have
!  affected the results in normal situations though.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  2/9/94 (D. Jahn)
!  Add 4th-order momentum advection.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  11/5/1995 (M. Xue)
!  Added vertical fourth order advection for w.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  9/1/98 (D. Weber and T. Chung)
!  Removed operators and merged loops to improve code efficiency.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mapfct   Map factors at scalar points
!
!    w        z component of velocity at a given time level (m/s)
!
!    ustr     u * rhostr
!    vstr     v * rhostr
!    wstr     w * rhostr
!
!  OUTPUT:
!
!    wforce   Acoustically inactive forcing terms in w-momentum
!             equation (kg/(m*s)**2).
!
!  WORK ARRAYS:
!
!    wadv     Advection term in w equation (kg/m**3)*(m/s)/s
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INCLUDE 'timelvls.inc'

  REAL :: mapfct(nx,ny)        ! Map factors at scalar points

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity at a given time level
                               ! (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity at a given time level
                               ! (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity at a given time level
                               ! (m/s)
  REAL :: ustr  (nx,ny,nz)     ! u*rhostr
  REAL :: vstr  (nx,ny,nz)     ! v*rhostr
  REAL :: wstr  (nx,ny,nz)     ! w*rhostr

  REAL :: wadv  (nx,ny,nz)     ! w-eqn advection term (kg/m**3)*(m/s)/s

  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: ibgn,iend,jbgn,jend
  REAL :: ueb, uwb, vsb, vnb
  REAL :: fourthirds,foursixth,onesix
  REAL :: tema

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'       ! Global model control parameters
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'bndry.inc'         ! Parameters for boundary conditions.
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  fourthirds = (4.0/3.0)
  foursixth = (4.0/6.0)
  onesix = (1.0/6.0)
!
!-----------------------------------------------------------------------
!
!  Calculate rho*u* dw/dx = avgx( avgz( ustr ) * difx( w ) )
!
!-----------------------------------------------------------------------

  tema = 0.25 * dxinv
  DO k=2,nz-1
    DO j=1,ny-1
      DO i=2,nx-2
        wadv(i,j,k)=tema*((ustr(i+1,j,k)+ustr(i+1,j,k-1))*              &
                          (w(i+1,j,k,2)-w(i,j,k,2))                     &
                         +(ustr(i,j,k)+ustr(i,j,k-1))*                  &
                          (w(i,j,k,2)-w(i-1,j,k,2)))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  At this point wadv contains avgx(avgz(u*rho)*dw/dx)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!    For 4th-order momentum advection, calculate the 4th-order
!    adjustment term for ustr*dw/dx:
!            avg2x( avgx(avgz(ustr)) * dif2x(w) )
!
!-----------------------------------------------------------------------

  IF (madvopt == 2 .OR. madvopt == 3) THEN

    tema = 0.125 * dxinv
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=2,nx-2
          tem2(i,j,k)=tema*((ustr(i,j,k-1)+ustr(i,j,k)) +               &
                            (ustr(i+1,j,k-1)+ustr(i+1,j,k))) *          &
                            (w(i+1,j,k,2)-w(i-1,j,k,2))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  At this point tem2 contains avgxz(ustr)*dw/dx
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!    In the case of periodic boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem2,nx,ny,nz,ebc,wbc,0,tem1)
    END IF

    CALL acct_interrupt(bc_acct)

    ibgn = 3
    iend = nx-3
    jbgn = 1
    jend = ny-1

    IF (ebc == 0) iend = nx-2

    IF (ebc == 2) THEN
      iend = nx-2
      IF (mp_opt == 0) THEN
        DO k=2,nz-1
        DO j=jbgn,jend
          tem2(nx-1,j,k) = tem2(2,j,k)
        END DO
        END DO
      END IF
    END IF

    IF (wbc == 0) ibgn = 2

    IF (wbc == 2) THEN
      ibgn = 2
      IF (mp_opt == 0) THEN
        DO k=2,nz-1
        DO j=jbgn,jend
          tem2(1,j,k) = tem2(nx-2,j,k)
        END DO
        END DO
      END IF
    END IF

!-----------------------------------------------------------------------
!
!    In the case of rigid boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (ebc == 1) THEN
      iend = nx-2
      DO k=2,nz-1
        DO j=jbgn,jend
          tem2(nx-1,j,k) = tem2(nx-2,j,k)
        END DO
      END DO
    END IF

    IF (wbc == 1) THEN
      ibgn = 2
      DO k=2,nz-1
        DO j=jbgn,jend
          tem2(1,j,k) = tem2(2,j,k)
        END DO
      END DO
    END IF

    CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Concatenate the 2nd- and 4th-order terms.
!
!-----------------------------------------------------------------------

    DO k=2,nz-1
      DO j=jbgn,jend
        DO i=ibgn,iend
          wadv(i,j,k) =  fourthirds*wadv(i,j,k)-                        &
                         onesix*(tem2(i+1,j,k)+tem2(i-1,j,k))
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Calculate rho*u* dw/dx on the east and west boundaries
!  using one-sided advection.
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=1,ny-1
      uwb=(ustr(1,j,k)+ustr(2,j,k)+ustr(1,j,k-1)+ustr(2,j,k-1))*0.25
      ueb=(ustr(nx-1,j,k)+ustr(nx,j,k)+ustr(nx-1,j,k-1)                 &
          +ustr(nx,j,k-1))*0.25
      IF (uwb < 0.0) THEN
        wadv(1,j,k)=uwb*(w(2,j,k,tpresent)-w(1,j,k,tpresent))*dxinv
      ELSE
        wadv(1,j,k)=rlxlbc*uwb*w(1,j,k,tpast)*dxinv
      END IF

      IF (ueb > 0.0) THEN
        wadv(nx-1,j,k)=ueb*                                             &
            (w(nx-1,j,k,tpresent)-w(nx-2,j,k,tpresent))*dxinv
      ELSE
        wadv(nx-1,j,k)=-rlxlbc*ueb*w(nx-1,j,k,tpast)*dxinv
      END IF
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate rho*v* dw/dy = avgy( avgz( vstr ) * dify ( w ) )
!
!-----------------------------------------------------------------------

  tema=0.25*dyinv
  DO k=2,nz-1
    DO j=2,ny-2
      DO i=1,nx-1
        tem2(i,j,k)=tema*((vstr(i,j+1,k)+vstr(i,j+1,k-1)) *             &
                          (w(i,j+1,k,2)-w(i,j,k,2))                     &
                         +(vstr(i,j,k)+vstr(i,j,k-1)) *                 &
                          (w(i,j,k,2)-w(i,j-1,k,2)))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  At this point tem2 contains avgz(v*rho)*dw/dy  second order term
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!    For 4th-order momentum advection, calculate the adjustment
!    term for vstr*dw/dy:
!            avg2y( avgz(avgy(vstr)) * dif2y(w) )
!
!-----------------------------------------------------------------------

  IF (madvopt == 2 .OR. madvopt == 3) THEN

    tema = 0.125*dyinv
    DO k=2,nz-1
      DO j=2,ny-2
        DO i=1,nx-1
          tem3(i,j,k)=tema*((vstr(i,j+1,k)+vstr(i,j,k)) +               &
                            (vstr(i,j+1,k-1)+vstr(i,j,k-1))) *          &
                            (w(i,j+1,k,2)-w(i,j-1,k,2))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  At this point tem3 contains avgyz(vstr)*dw/dy
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!    In the case of periodic boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dns(tem3,nx,ny,nz,nbc,sbc,0,tem1)
    END IF

    CALL acct_interrupt(bc_acct)

    ibgn = 1
    iend = nx-1
    jbgn = 3
    jend = ny-3

    IF (nbc == 0) jend = ny-2

    IF (nbc == 2) THEN
      jend = ny-2
      IF (mp_opt == 0) THEN
        DO k=2,nz-1
        DO i=ibgn,iend
          tem3(i,ny-1,k) = tem3(i,2,k)
        END DO
        END DO
      END IF
    END IF

    IF (sbc == 0) jbgn = 2

    IF (sbc == 2) THEN
      jbgn = 2
      IF (mp_opt == 0) THEN
        DO k=2,nz-1
        DO i=ibgn,iend
          tem3(i,1,k)=tem3(i,ny-2,k)
        END DO
        END DO
      END IF
    END IF

!-----------------------------------------------------------------------
!
!    In the case of rigid boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (nbc == 1) THEN
      jend = ny-2
      DO k=2,nz-1
        DO i=ibgn,iend
          tem3(i,ny-1,k) = tem3(i,ny-2,k)
        END DO
      END DO
    END IF

    IF (sbc == 1) THEN
      jbgn = 2
      DO k=2,nz-1
        DO i=ibgn,iend
          tem3(i,1,k)=tem3(i,2,k)
        END DO
      END DO
    END IF

    CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Concatenate the 2nd- and 4th- order terms.
!
!-----------------------------------------------------------------------

    DO k=2,nz-1
      DO j=jbgn,jend
        DO i=ibgn,iend
          tem2(i,j,k) =  fourthirds*tem2(i,j,k) -                       &
                         onesix*(tem3(i,j+1,k)+tem3(i,j-1,k))
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Calculate rho*v* dw/dy on the north and south boundaries
!  using one-sided advection.
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO i=1,nx-1
      vsb=(vstr(i,1,k)+vstr(i,2,k)+vstr(i,1,k-1)+vstr(i,2,k-1))*0.25
      vnb=(vstr(i,ny-1,k)+vstr(i,ny,k)+vstr(i,ny-1,k-1)                 &
          +vstr(i,ny,k-1))*0.25
      IF (vsb < 0.0) THEN
        tem2(i,1,k)=vsb*(w(i,2,k,tpresent)-w(i,1,k,tpresent))*dyinv
      ELSE
        tem2(i,1,k)=rlxlbc*vsb*w(i,1,k,tpast)*dyinv
      END IF

      IF (vnb > 0.0) THEN
        tem2(i,ny-1,k)=vnb*                                             &
            (w(i,ny-1,k,tpresent)-w(i,ny-2,k,tpresent))*dyinv
      ELSE
        tem2(i,ny-1,k)=-rlxlbc*vnb*w(i,ny-1,k,tpast)*dyinv
      END IF
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Add the w-equation horizontal advection terms, rho*u*dw/dx (wadv) and
!  rho*v*dw/dy (tem1).
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        wforce(i,j,k)=wforce(i,j,k) -                                   &
                       (wadv(i,j,k)+tem2(i,j,k)) * mapfct(i,j)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate rho*w* dw/dz = avgz( avgz( wstr ) * difz( w ) )
!
!-----------------------------------------------------------------------

  tema = 0.5*dzinv
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=tema*(wstr(i,j,k+1)+wstr(i,j,k))*                   &
                          (w(i,j,k+1,2)-w(i,j,k,2))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  tem2 contains avgz(w*rho)*dw/dz
!
!-----------------------------------------------------------------------

  IF (madvopt == 3) THEN

!-----------------------------------------------------------------------
!
!  If momentum advection 4th-order, calculate adjustment term:
!  avg2z(avg2z(wstr)*dif2z(w))
!
!-----------------------------------------------------------------------

    tema = 0.25*dzinv
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=tema*(wstr(i,j,k+1)+wstr(i,j,k-1))*               &
                            (w(i,j,k+1,2)-w(i,j,k-1,2))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Add the fourth order contribution to the 2nd order term and
!  add the w-equation vertical advection term rho*w*dw/dz (tem1)
!  to the horizontal advection terms.
!
!-----------------------------------------------------------------------

    DO k=3,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          wforce(i,j,k) = wforce(i,j,k) -                               &
                      foursixth*(tem1(i,j,k-1)+tem1(i,j,k)) +           &
                         onesix*(tem2(i,j,k+1)+tem2(i,j,k-1))
        END DO
      END DO
    END DO

    DO k=2,nz-1,nz-3
      DO j=1,ny-1
        DO i=1,nx-1
          wforce(i,j,k)=wforce(i,j,k) - (tem1(i,j,k)+tem1(i,j,k-1))*0.5
        END DO
      END DO
    END DO

  ELSE  ! perform second order advection....

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          wforce(i,j,k)=wforce(i,j,k) - (tem1(i,j,k)+tem1(i,j,k-1))*0.5
        END DO
      END DO
    END DO

  END IF


  RETURN
END SUBROUTINE advw
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADVPT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advpt(nx,ny,nz,ptprt,u,v,w,wcont,                            &
           rhostr,ptbar, mapfct,j3,j3inv,                               &
           ptadv,                                                       &
           ustr,vstr,wstr,tem1,tem2,tem3,tem4,mp_tem)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the advection of total potential temperature
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (M. Xue)
!  Added full documentation.
!
!  5/29/92 (K. Brewster)
!  Further facelift and order of argument list changed.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Added the terrain.
!
!  5/26/93 (M. Xue)
!  Corrected ptbar advection term and added w and j3 into the
!  argument list.
!
!  9/4/94 (D. Jahn)
!  Add tem4 to the call of advcts.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  3/8/96 (Ming Xue)
!  Fixed a bug in DO 510 i=1,nx-1, which was written as DO 510 i=i,nx-1.
!  Effectively, the vertical advection remained 2nd order with the error.
!
!  9/1/98 (D. Weber and T. Chung)
!  Removed operators and merged loops to improve code efficiency.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ptprt    Perturbation potential temperature at all time levels
!             (K)
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates at all time levels (m/s).
!    wcont    Contravariant vertical velocity (m/s)
!    rhostr   Base state air density (kg/m**3)
!    ptbar    Base state potential temperature (K)
!
!    mapfct   Map factors at scalar points
!
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    ptadv    Advection term of potential temperature eqn
!             (kg/m**3)*K/s
!
!  WORK ARRAYS:
!
!    ustr     u * rhostr
!    vstr     v * rhostr
!    wstr     w * rhostr
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INCLUDE 'timelvls.inc'

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)

  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)

  REAL :: mapfct(nx,ny)        ! Map factors at scalar points

  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               !  d(zp)/d(z)
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian
                               !  d(zp)/d(z)

  REAL :: ptadv(nx,ny,nz)      ! Potential temperature advection

  REAL :: ustr  (nx,ny,nz)     ! u * rhostr
  REAL :: vstr  (nx,ny,nz)     ! v * rhostr
  REAL :: wstr  (nx,ny,nz)     ! w * rhostr

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array

  REAL :: mp_tem(nx,ny,nz)     ! Temporary message passing array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,tlevel
  REAL :: tema
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'       ! Global model control parameters
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'bndry.inc'         ! Parameters for boundary conditions.
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tlevel=tpresent
!
!-----------------------------------------------------------------------
!
!  Calculate ustr=rhostr*u, vstr=rhostr*v, wstr=rhostr*w
!
!-----------------------------------------------------------------------

  CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        ustr(i,j,k)=u(i,j,k,tlevel)*tem1(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        vstr(i,j,k)=v(i,j,k,tlevel)*tem2(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx-1
        wstr(i,j,k)=wcont(i,j,k)*tem3(i,j,k)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Perturbation potential temperature advection
!
!-----------------------------------------------------------------------

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem3(i,j,k) = 0.0
      END DO
    END DO
  END DO

  CALL advcts(nx,ny,nz, ptprt,                                          &
              u(1,1,1,tlevel),v(1,1,1,tlevel),                          &
              ustr,vstr,wstr,tem3, mapfct,                              &
              ptadv, tem1,tem2,tem4,mp_tem)

  IF( ptsmlstp == 0 ) THEN

!-----------------------------------------------------------------------
!
!  Base state potential temperature advection.  This term is added to
!  the array ptadv to yield the total potential temperature advection.
!
!  ptbar is assumed to be independent of physical x and y, therefore
!  d(ptbar)/dx and d(ptbar)/dy for constant z are zero, the base state
!  advection is -w*d(ptbar)/dzp = -w/j3*d(ptbar)/dz.
!
!-----------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k) = rhostr(i,j,k)*j3inv(i,j,k)
        END DO
      END DO
    END DO

    tema=0.5*dzinv
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k)=tema*(tem2(i,j,k)+tem2(i,j,k-1))*w(i,j,k,2)*      &
                           (ptbar(i,j,k)-ptbar(i,j,k-1))
        END DO
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          ptadv(i,j,k)=ptadv(i,j,k)+0.5*(tem1(i,j,k+1)+tem1(i,j,k))
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE advpt
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ADVCTS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE advcts(nx,ny,nz,s,u,v,ustr,vstr,wstr,sbar, mapfct,           &
           sadv,                                                        &
           tem1,tem2,tem3,mp_tem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the scalar equation advection terms. These terms are
!  written in equivalent advection form.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (M. Xue)
!  Added full documentation.
!
!  4/9/93 (M. Xue & K. Brewster)
!  Some index bounds for operators and loop bounds were corrected.
!  Redundant calculations were done before at the boundaries that
!  caused out-of-bound calculations, which should not have
!  affected the results in normal situations though.
!
!  5/29/92 (K. Brewster)
!  Further facelift.
!
!  2/9/94 (D. Jahn)
!  Add 4th-order momentum advection.
!
!  9/1/94 (D. Weber & Y. Lu)
!  Cleaned up documentation
!
!  11/5/1995 (M. Xue)
!  Added vertical fourth order advection for scalars.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  3/8/96 (Ming Xue)
!  Fixed a bug in DO 510 i=1,nx-1, which was written as DO 510 i=i,nx-1.
!  Effectively, the vertical advection remained 2nd order with the error.
!
!  9/1/98 (D. Weber and T. Chung)
!  Removed operators and merged loops to improve code efficiency.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    mapfct   Map factors at scalar points
!
!    s        scalar field at a given time level (scalar units)
!
!    ustr     u * rhostr
!    vstr     v * rhostr
!    wstr     w * rhostr
!
!  OUTPUT:
!
!    sadv     Advection term in scalar equation (kg/m**3)*scalar
!             unit/s
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!
!    mp_tem   Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INCLUDE 'timelvls.inc'

  REAL :: mapfct(nx,ny)        ! Map factors at scalar points

  REAL :: s     (nx,ny,nz,nt)  ! a scalar field to be advected.
  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)

  REAL :: ustr  (nx,ny,nz)     ! u * rhostr
  REAL :: vstr  (nx,ny,nz)     ! v * rhostr
  REAL :: wstr  (nx,ny,nz)     ! w * rhostr
  REAL :: sbar  (nx,ny,nz)     ! A state of 's' towards which 's' is
                               ! relaxed on the inflow boundaries.

  REAL :: sadv  (nx,ny,nz)     ! advection term in scalar equation
                               ! (kg/m**3)*(scalar units)/s

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array

  REAL :: mp_tem(nx,ny,nz)     ! Temporary message passing array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: ibgn,iend,jbgn,jend
  REAL :: ueb, uwb, vsb, vnb
  REAL :: fourthirds,foursixth,onesix
  REAL :: tema
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Physical constants
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'bndry.inc'       ! Parameters for boundary conditions.
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  fourthirds = (4.0/3.0)
  foursixth = (4.0/6.0)
  onesix = (1.0/6.0)

!-----------------------------------------------------------------------
!
!  Calculate rho*u* ds/dx = avgx( ustr * difx ( s ) )
!
!-----------------------------------------------------------------------

  tema = dxinv*0.5
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=2,nx-2
        sadv(i,j,k)=tema*((s(i,j,k,2)-s(i-1,j,k,2))*ustr(i,j,k)         &
                         +(s(i+1,j,k,2)-s(i,j,k,2))*ustr(i+1,j,k))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!    For 4th-order scalar advection, calculate the 4th-order
!    adjustment term for ustr*ds/dx:
!            avg2x( avgx(ustr) * dif2x(s) )
!
!-----------------------------------------------------------------------

  IF (sadvopt == 2.OR.sadvopt == 3) THEN

    tema = 0.25*dxinv
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=2,nx-2
          tem2(i,j,k)=tema*(ustr(i+1,j,k)+ustr(i,j,k))*                 &
                           (s(i+1,j,k,2)-s(i-1,j,k,2))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!    In the case of periodic boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem2,nx,ny,nz,ebc,wbc,0,tem1)
    END IF

    CALL acct_interrupt(bc_acct)

    ibgn = 3
    iend = nx-3
    jbgn = 1
    jend = ny-1

    IF (ebc == 0) iend = nx-2

    IF (ebc == 2) THEN
      iend = nx-2
      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO j=jbgn,jend
          tem2(nx-1,j,k) = tem2(2,j,k)
        END DO
        END DO
      END IF
    END IF

    IF (wbc == 0) ibgn = 2

    IF (wbc == 2) THEN
      ibgn = 2
      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO j=jbgn,jend
          tem2(1,j,k) = tem2(nx-2,j,k)
        END DO
        END DO
      END IF
    END IF

!-----------------------------------------------------------------------
!
!    In the case of rigid boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (ebc == 1) THEN
      iend = nx-2
      DO k=2,nz-2
        DO j=jbgn,jend
          tem2(nx-1,j,k) = tem2(nx-2,j,k)
        END DO
      END DO
    END IF

    IF (wbc == 1) THEN
      ibgn = 2
      DO k=2,nz-2
        DO j=jbgn,jend
          tem2(1,j,k) = tem2(2,j,k)
        END DO
      END DO
    END IF

    CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Concatenate the 2nd- and 4th- order terms for ustr*ds/dx
!
!-----------------------------------------------------------------------

    DO  k=2,nz-2
      DO  j=jbgn,jend
        DO  i=ibgn,iend
          sadv(i,j,k) = fourthirds*sadv(i,j,k) -                        &
                        onesix*(tem2(i+1,j,k)+tem2(i-1,j,k))
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Calculate rho*u* ds/dx on the east and west boundaries
!  using one-sided advection.
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=1,ny-1
      uwb=(ustr(1,j,k)+ustr(2,j,k))*0.5
      ueb=(ustr(nx-1,j,k)+ustr(nx,j,k))*0.5
      IF (uwb < 0.0) THEN
        sadv(1,j,k)=uwb*(s(2,j,k,tpresent)-s(1,j,k,tpresent))*dxinv
      ELSE
        sadv(1,j,k)=rlxlbc*uwb*(s(1,j,k,tpast)-sbar(1,j,k))*dxinv
      END IF

      IF (ueb > 0.0) THEN
        sadv(nx-1,j,k)=ueb*                                             &
            (s(nx-1,j,k,tpresent)-s(nx-2,j,k,tpresent))*dxinv
      ELSE
        sadv(nx-1,j,k)=-rlxlbc*ueb*                                     &
            (s(nx-1,j,k,tpast)-sbar(nx-1,j,k))*dxinv
      END IF
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate rho*v* ds/dy = avgy( vstr * dify( s ) )
!
!-----------------------------------------------------------------------

  tema = 0.5*dyinv
  DO k=2,nz-2
    DO j=2,ny-2
      DO i=1,nx-1
        tem1(i,j,k)=tema*((s(i,j,k,2)-s(i,j-1,k,2))*vstr(i,j,k)         &
                         +(s(i,j+1,k,2)-s(i,j,k,2))*vstr(i,j+1,k))
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!    For 4th-order scalar advection, calculate the 4th-order
!    adjustment term for vstr*ds/dy:
!            avg2y( avgy(vstr) * dif2y(s) )
!
!-----------------------------------------------------------------------

  IF (sadvopt == 2.OR.sadvopt == 3) THEN

    tema = 0.25*dyinv
    DO k=2,nz-2
      DO j=2,ny-2
        DO i=1,nx-1
          tem3(i,j,k)=tema*(vstr(i,j+1,k)+vstr(i,j,k))*                 &
                            (s(i,j+1,k,2)-s(i,j-1,k,2))
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!    In the case of periodic boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dns(tem3,nx,ny,nz,nbc,sbc,0,mp_tem)
    END IF

    CALL acct_interrupt(bc_acct)

    ibgn = 1
    iend = nx-1
    jbgn = 3
    jend = ny-3

    IF (nbc == 0) jend = ny-2

    IF (nbc == 2) THEN
      jend = ny-2
      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,ny-1,k) = tem3(i,2,k)
        END DO
        END DO
      END IF
    END IF

    IF (sbc == 0) jbgn = 2

    IF (sbc == 2) THEN
      jbgn = 2
      IF (mp_opt == 0) THEN
        DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,1,k) = tem3(i,ny-2,k)
        END DO
        END DO
      END IF
    END IF

!-----------------------------------------------------------------------
!
!    In the case of rigid boundary condition, need to calculate
!    the 4th-order term for one additional point near boundary.
!
!-----------------------------------------------------------------------

    IF (nbc == 1) THEN
      jend = ny-2
      DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,ny-1,k) = tem3(i,ny-2,k)
        END DO
      END DO
    END IF

    IF (sbc == 1) THEN
      jbgn = 2
      DO k=2,nz-2
        DO i=ibgn,iend
          tem3(i,1,k) = tem3(i,2,k)
        END DO
      END DO
    END IF

    CALL acct_stop_inter

!-----------------------------------------------------------------------
!
!  Concatenate the 2nd- and 4th-order terms for vstr*ds/dy.
!
!-----------------------------------------------------------------------

    DO k=2,nz-2
      DO j=jbgn,jend
        DO i=ibgn,iend
          tem1(i,j,k) =  fourthirds*tem1(i,j,k) -                       &
                         onesix*(tem3(i,j+1,k)+tem3(i,j-1,k))
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
!  Calculate rho*v* ds/dy on the north and south boundaries
!  using one-sided advection.
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO i=1,nx-1
      vsb=(vstr(i,1,k)+vstr(i,2,k))*0.5
      vnb=(vstr(i,ny-1,k)+vstr(i,ny,k))*0.5
      IF (vsb < 0.0) THEN
        tem1(i,1,k)=vsb*(s(i,2,k,tpresent)-s(i,1,k,tpresent))*dyinv
      ELSE
        tem1(i,1,k)=rlxlbc*vsb*(s(i,1,k,tpast)-sbar(i,1,k))*dyinv
      END IF

      IF (vnb > 0.0) THEN
        tem1(i,ny-1,k)=vnb*                                             &
            (s(i,ny-1,k,tpresent)-s(i,ny-2,k,tpresent))*dyinv
      ELSE
        tem1(i,ny-1,k)=-rlxlbc*vnb*                                     &
            (s(i,ny-1,k,tpast)-sbar(i,ny-1,k))*dyinv
      END IF
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Add the scalar equation horizontal advection terms, rho*u*ds/dx (sadv)
!  and  rho*v*ds/dy (tem1).
!
!-----------------------------------------------------------------------

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        sadv(i,j,k)=(sadv(i,j,k)+tem1(i,j,k)) * mapfct(i,j)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Calculate rho*w* ds/dz = avgz( wstr * difz( s ) )
!
!-----------------------------------------------------------------------

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k)=dzinv*(s(i,j,k,2)-s(i,j,k-1,2))*wstr(i,j,k)
      END DO
    END DO
  END DO

  IF (sadvopt == 3) THEN

!-----------------------------------------------------------------------
!
!  For 4th-order scalar advection, calculate the 4th-order
!  adjustment term for wstr*ds/dz: avg2z( avgz(wstr) * dif2z(s) )
!  and combine it with the second order term.
!
!  Note that we apply the 4th order scheme only upto the second point
!  from the top and bottom boundary. We do not give special treatment
!  to the periodic or wall boundary cases as we do for the side
!  boundaries.
!
!  Add the scalar equation vertical advection term rho*w*ds/dz (tem2)
!  to the horizontal advection terms.
!
!-----------------------------------------------------------------------

    tema = 0.25*dzinv
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          tem3(i,j,k)=tema*(wstr(i,j,k+1)+wstr(i,j,k))*                 &
                            (s(i,j,k+1,2)-s(i,j,k-1,2))
        END DO
      END DO
    END DO

    DO k=3,nz-3
      DO j=1,ny-1
        DO i=1,nx-1
          sadv(i,j,k) = sadv(i,j,k)+                                    &
                      foursixth*(tem1(i,j,k+1)+tem1(i,j,k)) -           &
                         onesix*(tem3(i,j,k+1)+tem3(i,j,k-1))
        END DO
      END DO
    END DO

    DO k=2,nz-2,nz-4
      DO j=1,ny-1
        DO i=1,nx-1
          sadv(i,j,k)= sadv(i,j,k) + (tem1(i,j,k+1)+tem1(i,j,k))*0.5
        END DO
      END DO
    END DO

  ELSE     !  perform second order advection...

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          sadv(i,j,k)= sadv(i,j,k) + (tem1(i,j,k+1)+tem1(i,j,k))*0.5
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE advcts
