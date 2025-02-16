!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE NESTBDT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE nestbdt(nx,ny,nz, a,                                         &
           iwb,ieb,jsb,jnb,kbb,ktb, dtbig,                              &
           adteb,adtwb,adtnb,adtsb)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the time tendency of array (a) on the lateral boundaries
!  for a nested grid.
!
!  abdt = (a(tpresent)-a(tpast))/dtbig
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/10/92.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east /west )
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    a        Mixing ratio for one of the water variables (kg/kg)
!    iwb,ieb,jsb,jnb,kbb,ktb
!             Indices specifying domain boundaries of nested grid
!    dtbig    Big time step size (s)
!
!  OUTPUT:
!
!    adteb    Time tendency of a at east  boundary (kg/(kg*s))
!    adtwb    Time tendency of a at west  boundary (kg/(kg*s))
!    adtnb    Time tendency of a at north boundary (kg/(kg*s))
!    adtsb    Time tendency of a at south boundary (kg/(kg*s))
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

  INTEGER :: nx,ny,nz     ! The number of grid points in 3 directions

  REAL :: dtbig           ! The big time step size (s)

  REAL :: a(nx,ny,nz,nt)  ! Mixing ratio for one of the water
                          ! variables (kg/kg)
  INTEGER :: iwb,ieb,jsb,jnb,kbb,ktb
                           ! Index for the variable boundaries.

  REAL :: adteb (ny,nz)   ! T-tendency of a at e-boundary (kg/(kg*s))
  REAL :: adtwb (ny,nz)   ! T-tendency of a at w-boundary (kg/(kg*s))
  REAL :: adtnb (nx,nz)   ! T-tendency of a at n-boundary (kg/(kg*s))
  REAL :: adtsb (nx,nz)   ! T-tendency of a at s-boundary (kg/(kg*s))
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL :: rdtbig
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  rdtbig = 1.0/dtbig

  DO k=kbb,ktb
    DO j=jsb,jnb
      adtwb(j,k)=(a(iwb,j,k,tpresent)-a(iwb,j,k,tpast))*rdtbig
      adteb(j,k)=(a(ieb,j,k,tpresent)-a(ieb,j,k,tpast))*rdtbig
    END DO
  END DO
!
  DO k=kbb,ktb
    DO i=iwb,ieb
      adtsb(i,k)=(a(i,jsb,k,tpresent)-a(i,jsb,k,tpast))*rdtbig
      adtnb(i,k)=(a(i,jnb,k,tpresent)-a(i,jnb,k,tpast))*rdtbig
    END DO
  END DO
!
  RETURN
END SUBROUTINE nestbdt
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BKWSMLDT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bkwsmldt(nx,ny,nz, u,v,ubar,vbar,                            &
           udteb,udtwb,vdtnb,vdtsb)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the time tendences of u and v on the lateral boundaries
!  using Klemp and Wilhelmson type radiation boundary condition
!  implemented inside the small time steps.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  12/28/1994.
!
!  MODIFICATION HISTORY
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east /west )
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x-component of velocity at future time levels (m/s)
!    v        y-component of velocity at future time levels (m/s)
!    ubar     Base state u-velocity (m/s)
!    vbar     Base state v-velocity (m/s)
!
!  OUTPUT:
!
!    udteb    Time tendency of u at east  boundary (m/s**2)
!    udtwb    Time tendency of u at west  boundary (m/s**2)
!    vdtnb    Time tendency of v at north boundary (m/s**2)
!    vdtsb    Time tendency of v at south boundary (m/s**2)
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

  INTEGER :: nx,ny,nz     ! The number of grid points in 3 directions

  REAL :: u   (nx,ny,nz)  ! U-velocity at time=tfuture (m/s).
  REAL :: v   (nx,ny,nz)  ! Total v-velocity (m/s)
  REAL :: ubar(nx,ny,nz)  ! Base state U-velocity (m/s).
  REAL :: vbar(nx,ny,nz)  ! Base state v-velocity (m/s)

  REAL :: udteb (ny,nz)   ! Time-tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)   ! Time-tendency of u at w-boundary (m/s**2)
  REAL :: vdtnb (nx,nz)   ! Time-tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)   ! Time-tendency of v at s-boundary (m/s**2)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL :: cphase
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'globcst.inc'     ! Global model control parameters
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
!  Calculate the time tendency of u on the west boundary (udtwb)
!
!-----------------------------------------------------------------------
!
  IF(wbc == 4) THEN             ! Radiation boundary condition.

    DO k=1,nz-1
      DO j=1,ny-1

        cphase = u(1,j,k)-c_phase

        udtwb(j,k)=(-MIN(cphase,0.0)*(u(2,j,k)-u(1,j,k))                &
                   -rlxlbc*MAX(cphase,0.0)*(u(1,j,k)-ubar(1,j,k))       &
                   )*dxinv

      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the time tendency of u on the east boundary (udteb)
!
!-----------------------------------------------------------------------
!
  IF(ebc == 4) THEN             ! Radiation boundary condition.

    DO k=1,nz-1
      DO j=1,ny-1

        cphase = u(nx,j,k)+c_phase

        udteb(j,k)=(-MAX(cphase,0.0)*(u(nx,j,k)-u(nx-1,j,k))            &
                   +rlxlbc*MIN(cphase,0.0)*(u(nx,j,k)-ubar(nx,j,k))     &
                   )*dxinv

      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the time tendency of v on the south boundary (vdtsb)
!
!-----------------------------------------------------------------------
!
  IF(sbc == 4) THEN             ! Radiation boundary condition.

    DO k=1,nz-1
      DO i=1,nx-1

        cphase = v(i,1,k)-c_phase

        vdtsb(i,k)=(-MIN(cphase,0.0)*(v(i,2,k)-v(i,1,k))                &
                   -rlxlbc*MAX(cphase,0.0)*(v(i,1,k)-vbar(i,1,k))       &
                   )*dyinv

      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the time tendency of v on the north boundary (vdtnb)
!
!-----------------------------------------------------------------------
!
  IF(nbc == 4) THEN             ! Radiation boundary condition.

    DO k=1,nz-1
      DO i=1,nx-1

        cphase = v(i,ny,k)+c_phase

        vdtnb(i,k)=(-MAX(cphase,0.0)*(v(i,ny,k)-v(i,ny-1,k))            &
                   +rlxlbc*MIN(cphase,0.0)*(v(i,ny,k)-vbar(i,ny,k))     &
                   )*dyinv

      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE bkwsmldt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BDTU                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bdtu(nx,ny,nz,dtbig1, u,ubar,udteb,udtwb)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the time tendency of u on the lateral boundaries
!  using the radiation boundary condition of Klemp-Wilhemlson or
!  Klemp-Lilly (1980)/Durran (1983).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber and Ming Xue.
!  3/27/95.
!
!  Includes vertical averaging of the computed orlanski phase
!  speeds following Klemp-Lilly (1978)/Durran (1983) approach
!  (rbcopt=3). Also includes the Klemp-Wilhelmson condition performed
!  on the big time step(rbcopt=2.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east /west )
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x-component of velocity at all time levels (m/s).
!
!  OUTPUT:
!
!    udteb    Time tendency of u at east  boundary (m/s**2)
!    udtwb    Time tendency of u at west  boundary (m/s**2)
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

  INTEGER :: nx,ny,nz     ! The number of grid points in 3 directions

  REAL :: dtbig1          ! The big time step size (s)

  REAL :: u (nx,ny,nz,nt) ! Total u-velocity (m/s).
  REAL :: ubar(nx,ny,nz)

  REAL :: udteb (ny,nz)  ! T-tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)  ! T-tendency of u at w-boundary (m/s**2)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: j, k
  REAL :: eps, tem
  REAL :: avguwest,avgueast

  REAL :: cdtdx           ! Estimated phase speed times dt/dx .
  REAL :: rdtbig          ! rdtbig = 1/dtbig1
  REAL :: cflimit         ! added for durran condition
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'globcst.inc'     ! Global model control parameters
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  eps=1.0E-30
  rdtbig = 1.0/dtbig1
  cflimit = 1.0
!
!-----------------------------------------------------------------------
!
!  Calculate the time tendency of u on the west boundary (udtwb)
!  The variable eps is added to the denominator to avoid division by
!  zero.
!
!  For the west boundary, cdtdx must be negative but not less than -1
!
!-----------------------------------------------------------------------
!
  IF(wbc == 4) THEN             ! Radiation boundary condition.

    DO j=1,ny-1
      avguwest = 0.0
      DO k=2,nz-2

        IF( rbcopt == 2) THEN     ! Klemp-Wilhelmson condition

          udtwb(j,k)=( -MIN(u(1,j,k,3)-c_phase, 0.0)*                   &
                     (u(2,j,k,3)-u(1,j,k,3))                            &
                     -rlxlbc*MAX(u(1,j,k,3)-c_phase,0.0)                &
                     *(u(1,j,k,2)-ubar(1,j,k)) )*dxinv

        ELSE IF( rbcopt == 3.OR.rbcopt == 4) THEN ! Orlanski Method
                                                  ! to compute cdtdx...

          tem  =-u(2,j,k,3)-u(2,j,k,1)+2.*u(3,j,k,2)
          cdtdx=(u(2,j,k,1)-u(2,j,k,3))/(SIGN(eps,tem)+tem)
          cdtdx=MAX(-cflimit, MIN(cdtdx,cflimit) )
          IF(cdtdx < 0.0) avguwest = avguwest + cdtdx

        END IF

        IF(rbcopt == 3) THEN     ! Complete the Orlanski computation
          udtwb(j,k)=( -MIN(cdtdx,0.0)*(u(2,j,k,3)-u(1,j,k,2))/         &
                     (1.-MIN(0.0,cdtdx))                                &
                  -rlxlbc*MAX(cdtdx,0.0)*(u(1,j,k,2)-ubar(1,j,k))       &
                  )*rdtbig
        END IF

      END DO

      IF( rbcopt == 4)THEN ! Apply Klemp-Lilly/Durran condition...

        cdtdx = avguwest/(nz-3)

        DO k=2,nz-2
          udtwb(j,k)=( -MIN(cdtdx,0.0)*(u(2,j,k,3)-u(1,j,k,2))/         &
                     (1.-MIN(0.0,cdtdx))                                &
                     -rlxlbc*MAX(cdtdx,0.0)*(u(1,j,k,2)-ubar(1,j,k))    &
                     )*rdtbig
        END DO

      END IF

    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the time tendency of u on the east boundary (udteb)
!  The variable eps is added to the denominator
!  to avoid division by zero.
!
!  For the east boundary, cdtdx must be positive but not larger
!  than 1
!
!-----------------------------------------------------------------------
!
  IF(ebc == 4) THEN             ! Radiation boundary condition.

    DO j=1,ny-1
      avgueast = 0.0
      DO k=2,nz-2

        IF( rbcopt == 2) THEN     ! Klemp-Wilhelmson condition

          udteb(j,k)=(-MAX(u(nx,j,k,3)+c_phase, 0.0)*                   &
                     (u(nx,j,k,3)-u(nx-1,j,k,3))                        &
                     +rlxlbc*MIN(u(nx,j,k,3)+c_phase,0.)                &
                     *(u(nx,j,k,2)-ubar(nx,j,k)) )*dxinv

        ELSE IF( rbcopt == 3.OR.rbcopt == 4) THEN ! Orlanski Method
                                                  ! to compute cdtdx...

          tem  =u(nx-1,j,k,3)+u(nx-1,j,k,1)-2.*u(nx-2,j,k,2)
          cdtdx=(u(nx-1,j,k,1)-u(nx-1,j,k,3))/(SIGN(eps,tem)+tem)
          cdtdx=MAX(-cflimit, MIN(cdtdx,cflimit))
          IF(cdtdx > 0.0) avgueast = avgueast + cdtdx

        END IF

        IF(rbcopt == 3) THEN     ! Complete the Orlanski computation
          udteb(j,k)=( MAX(cdtdx,0.0)*(u(nx-1,j,k,3)-u(nx,j,k,2))/      &
                   (1.+MAX(0.0,cdtdx))                                  &
                   +rlxlbc*MIN(cdtdx,0.0)*(u(nx,j,k,2)-ubar(nx,j,k))    &
                   )*rdtbig

        END IF

      END DO

      IF( rbcopt == 4)THEN ! Apply Klemp-Lilly/Durran condition...

        cdtdx = avgueast/(nz-3)

        DO k=2,nz-2
          udteb(j,k)=( MAX(cdtdx,0.0)*(u(nx-1,j,k,3)-u(nx,j,k,2))/      &
                   (1.+MAX(0.0,cdtdx))                                  &
                   +rlxlbc*MIN(cdtdx,0.0)*(u(nx,j,k,2)-ubar(nx,j,k))    &
                   )*rdtbig
        END DO

      END IF

    END DO

  END IF

  RETURN
END SUBROUTINE bdtu

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BDTV                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bdtv(nx,ny,nz,dtbig1, v,vbar,vdtnb,vdtsb)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the time tendency of v on the lateral boundaries
!  using the radiation boundary condition of Klemp-Wilhemlson or
!  Klemp-Lilly (1980)/Durran (1983).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber and Ming Xue.
!  3/27/95.
!
!  Includes vertical averaging of the computed orlanski phase
!  speeds following Klemp-Lilly (1978)/Durran (1983) approach
!  (rbcopt=3). Also includes the Klemp-Wilhelmson condition performed
!  on the big time step(rbcopt=2.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east /west )
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    v        y-component of velocity at all time levels (m/s).
!
!  OUTPUT:
!
!    vdtnb    Time tendency of v at north boundary (m/s**2)
!    vdtsb    Time tendency of v at south boundary (m/s**2)
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

  INTEGER :: nx,ny,nz     ! The number of grid points in 3 directions

  REAL :: dtbig1          ! The big time step size (s)

  REAL :: v (nx,ny,nz,nt) ! Total v-velocity (m/s).
  REAL :: vbar(nx,ny,nz)

  REAL :: vdtnb (nx,nz)   ! T-tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)   ! T-tendency of v at s-boundary (m/s**2)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, k
  REAL :: eps, tem
  REAL :: avgvsouth,avgvnorth

  REAL :: cdtdy           ! Estimated phase speed times dt/dy .
  REAL :: rdtbig          ! rdtbig = 1/dtbig1
  REAL :: cflimit         ! added for durran condition...
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'globcst.inc'     ! Global model control parameters
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  eps=1.0E-30
  rdtbig = 1.0/dtbig1
  cflimit = 1.0
!
!-----------------------------------------------------------------------
!
!  Calculate the time tendency of v on the south boundary (vdtsb)
!  The variable eps is added to the denominator to avoid division by
!  zero.
!
!  For the south boundary, cdtdy must be negative but not less
!  than -1
!
!-----------------------------------------------------------------------
!
  IF(sbc == 4) THEN             ! Radiation boundary condition.

    DO i=1,nx-1
      avgvsouth = 0.0
      DO k=2,nz-2

        IF( rbcopt == 2) THEN     ! Klemp-Wilhelmson condition

          vdtsb(i,k)=(-MIN(v(i,1,k,3)-c_phase, 0.0)*                    &
                     (v(i,2,k,3)-v(i,1,k,3))                            &
                     -rlxlbc*MAX(v(i,1,k,3)-c_phase,0.)                 &
                     *(v(i,1,k,2)-vbar(i,1,k)) )*dyinv

        ELSE IF( rbcopt == 3.OR.rbcopt == 4) THEN ! Orlanski Method
                                                  ! to compute cdtdx...

          tem  =-v(i,2,k,3)-v(i,2,k,1)+2.*v(i,3,k,2)
          cdtdy=(v(i,2,k,1)-v(i,2,k,3))/(SIGN(eps,tem)+tem)
          cdtdy=MAX(-cflimit, MIN(cdtdy,cflimit))
          IF(cdtdy < 0.0) avgvsouth = avgvsouth + cdtdy

        END IF

        IF(rbcopt == 3) THEN     ! Complete the Orlanski computation
          vdtsb(i,k)=( MIN(cdtdy,0.0)*(-v(i,2,k,3)+v(i,1,k,2))/         &
                     (1.-MIN(0.0,cdtdy))                                &
                     -rlxlbc*MAX(cdtdy,0.0)                             &
                     *(v(i,1,k,2)-vbar(i,1,k)))*rdtbig
        END IF

      END DO

      IF( rbcopt == 4)THEN ! Apply Klemp-Lilly/Durran condition...

        cdtdy = avgvsouth/(nz-3)

        DO k=2,nz-2
          vdtsb(i,k)=( MIN(cdtdy,0.0)*(-v(i,2,k,3)+v(i,1,k,2))/         &
                     (1.-MIN(0.0,cdtdy))                                &
                     -rlxlbc*MAX(cdtdy,0.0)                             &
                     *(v(i,1,k,2)-vbar(i,1,k)))*rdtbig
        END DO

      END IF

    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the time tendency of v on the north boundary (vdtnb)
!  The variable eps is added to the denominator
!  to avoid division by zero.
!
!  For the north boundary, cdtdy must be positive but not larger
!  than 1
!
!-----------------------------------------------------------------------
!
  IF(nbc == 4) THEN             ! Radiation boundary condition.

    DO i=1,nx-1
      avgvnorth = 0.0
      DO k=2,nz-2

        IF( rbcopt == 2) THEN     ! Klemp-Wilhelmson condition

          vdtnb(i,k)=(-MAX(v(i,ny,k,3)+c_phase, 0.0)*                   &
                     (v(i,ny,k,3)-v(i,ny-1,k,3))                        &
                     +rlxlbc*MIN(v(i,ny,k,3)+c_phase,0.)                &
                     *(v(i,ny,k,2)-vbar(i,ny,k)) )*dyinv


        ELSE IF( rbcopt == 3.OR.rbcopt == 4) THEN ! Orlanski Method
                                                  ! to compute cdtdx...

          tem  =v(i,ny-1,k,3)+v(i,ny-1,k,1)-2.*v(i,ny-2,k,2)
          cdtdy=(v(i,ny-1,k,1)-v(i,ny-1,k,3))/(SIGN(eps,tem)+tem)
          cdtdy=MAX(-cflimit, MIN(cdtdy,cflimit))
          IF(cdtdy > 0.0) avgvnorth = avgvnorth + cdtdy

        END IF

        IF(rbcopt == 3) THEN     ! Complete the Orlanski computation
          vdtnb(i,k)=( MAX(cdtdy,0.0)*(v(i,ny-1,k,3)-v(i,ny,k,2))/      &
                   (1.+MAX(0.0,cdtdy))                                  &
                   +rlxlbc*MIN(cdtdy,0.0)*(v(i,ny,k,2)-vbar(i,ny,k))    &
                   )*rdtbig
        END IF

      END DO

      IF( rbcopt == 4)THEN ! Apply Klemp-Lilly/Durran condition...

        cdtdy = avgvnorth/(nz-3)

        DO k=2,nz-2
          vdtnb(i,k)=( MAX(cdtdy,0.0)*(v(i,ny-1,k,3)-v(i,ny,k,2))/      &
                   (1.+MAX(0.0,cdtdy))                                  &
                   +rlxlbc*MIN(cdtdy,0.0)*(v(i,ny,k,2)-vbar(i,ny,k))    &
                   )*rdtbig
        END DO

      END IF

    END DO

  END IF

  RETURN
END SUBROUTINE bdtv
