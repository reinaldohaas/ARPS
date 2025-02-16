!
!##################################################################
!##################################################################
!######                                                      ######
!######         DUMMY SUBROUTINES FOR ASSIMILATION           ######
!######                                                      ######
!######                Copyright (c) 1993                    ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE assimcon(nx,ny,nz,x,y,z,zp,                             &
                  ubar,vbar,pbar,ptbar,rhostr,qvbar,               &
                  u,v,w,ptprt,pprt,qv,qscalar,                     &
                  j1,j2,j3,                                        &
                  tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,    &
                  tem10,tem11)

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE         ! Force explicit declarations

  INCLUDE 'timelvls.inc'

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point(m).
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point(m).
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid(m).
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid(m).
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nt)  ! should be (nx,ny,ny,nt,nscalar)

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/d(x)
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/d(y)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! d(zp)/d(z)

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

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(a,a/a/a)')                                                  &
      'Dummy subroutine ASSIMCON was called. ',                         &
      'Check input parameter assimopt. ',                               &
      'If assimopt=1, the real subroutine should be linked.',           &
      'Program stopped in dummy ASSIMCON in file assimdummy.f'

  CALL arpsstop('arpsstop called from assimcon',1)

END SUBROUTINE assimcon
!
!##################################################################
!##################################################################
!######                                                      ######
!######         DUMMY SUBROUTINES FOR ASSIMDRIVi             ######
!######                                                      ######
!######                Copyright (c) 1993                    ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE assimdriv(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                &
                   u,v,w,wcont,ptprt,pprt,qv,qscalar,                &
                   tke,kmh,kmv,mapfct,                               &
                   ubar,vbar,ptbar,pbar,rhostr,qvbar,                &
                   j1,j2,j3,j3soil,                                  &
                   tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,     &
                   tem10,tem11,tem12,                                &
                   tem13,tem14,tem15,tem16)

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INCLUDE 'timelvls.inc'

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels

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

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nt)  ! should be (nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent kinetic energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/d(x)
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/d(y)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! d(zp)/d(z)
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian
                               ! d(zpsoil)/d(z)

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
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(a,a/a/a)')                                                  &
      'Dummy subroutine ASSIMDRIV was called. ',                        &
      'Check input parameter assimopt. ',                               &
      'If assimopt=1, the real subroutine should be linked.',           &
      'Program stopped in dummy ASSIMDRIV in file assimdummy.f'

  CALL arpsstop('arpsstop called from assimdriv ',1)

END SUBROUTINE assimdriv
!
!##################################################################
!##################################################################
!######                                                      ######
!######         DUMMY SUBROUTINES FOR ASSIMPTPR              ######
!######                                                      ######
!######                Copyright (c) 1993                    ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE  assimptpr(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,             &
                   u,v,w,ptprt,pprt,qv,qscalar,                    &
                   tke,kmh,kmv,                                    &
                   ubar,vbar,ptbar,pbar,qvbar,rhostr,              &
                   uforce,vforce,wforce,                           &
                   j1,j2,j3,j3soil)

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INCLUDE 'timelvls.inc'

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels

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

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nt)  ! should be (nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent kinetic energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)
                               ! uforce= -uadv + umix + ucorio
  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)
                               ! vforce= -vadv + vmix + vcorio
  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
                               ! wforce= -wadv + wmix + wbuoy

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/d(x)
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/d(y)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! d(zp)/d(z)
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian
                               ! d(zpsoil)/d(z)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  WRITE (6,'(a,a,a)') 'You are calling a dummy function. ',             &
                    'The model stopped here in assimdummy.f.'
  CALL arpsstop('arpsstop called from assimdriv',1)
END SUBROUTINE assimptpr
