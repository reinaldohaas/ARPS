!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INIBASE                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE inibase(nx,ny,nz,                                            &
           ubar,vbar,ptbar,pbar,ptbari,pbari,                           &
           rhostr,rhostri,qvbar,                                        &
           x,y,z,zp,j3, rhobar, zs,zuv,pibar, tem1, tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the base-state variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  10/05/92 (M. Xue)
!  zero gradient in ptbar at the bottom and top boundaries are
!  no longer imposed.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  10/31/94 (E. Adlerman)
!  Added option for dewpoint input on sounding and for
!  inibasopt = 5, i.e. Weisman&Klemp analytical thermodynamic
!  profile.
!
!  12/22/94 (Yuhe Liu)
!  Added option for base state wind input for inibasopt != 1 and
!  viniopt = 1 and 2.
!
!  5/10/95 (M. Xue)
!  Corrected a bug with the pibar calculation for inibasopt=4 case.
!
!  5/10/95 (M. Xue)
!  Changed lower boundary condition for ptbar and qvbar to constant
!  gradient extrapolation. rhobar, pbar and pibar at k=1 are obtained
!  from hydrostatic equation.
!
!  7/6/95 (M. Xue)
!  Zero gradient condition restored for ptbar and qvbar at the
!  bottom boundary. The extrapolated condition was introducing
!  undesirable turbulent heat and moisture fluxes at the ground.
!
!  2/11/1997 (M. Xue and D. Weber)
!  Corrected a problem with base-state initialization option 4,
!  the case of constant static stability.
!
!  8/4/1997 (M. Xue and D. Weber)
!  Fixed a problem with pibar calculation for option 4, introduced
!  on 2/11/1997.
!
!  5/1/1998 (G. Bassett and D. Weber)
!  Moved sounding output to subroutine writesnd.
!
!  9/15/1998 (D. Weber)
!  Added ptbari, pbari, and rhostri (inverse quantities).
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse base state potential temperature (K)
!    pbari    Inverse base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    pibar    Base state
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!  WORK ARRAYS:
!
!    rhobar   Temporary work array.
!    zs       Temporary work array.
!    zuv      Temporary work array.
!    tem1     Work array
!    tem2     Work array
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

  INCLUDE 'mp.inc'          ! Message passing parameters.


  INTEGER :: nx,ny,nz          ! The number of grid points in 3
                               ! directions

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water specific humidity
                               ! (kg/kg).
  REAL :: pibar (nx,ny,nz)     ! Base state Exner function
  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of staggered grid.
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! as d(zp)/dz.

  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: zs    (nx,ny,nz)     ! Physical coordinate height at scalar
                               ! point (m)
  REAL :: zuv   (nx,ny,nz)     ! Physical coordinate height at u or v
                               ! point (m)

  REAL :: tem1  (nx,ny,nz)
  REAL :: tem2  (nx,ny,nz)

  INTEGER :: lvlprof
  REAL :: depthp

  PARAMETER( lvlprof=601, depthp = 3.0E4 )

!
!-----------------------------------------------------------------------
!
!  lvlprof:
!
!  The ARPS interpolates unevenly-spaced data from a vertical
!  sounding to evenly-spaced data for use in defining the base
!  state atmosphere.  In this process, an intermediate sounding is
!  generated at evenly-spaced altitudes, with the accuracy of the
!  associated interpolation controlled by the parameter lvlprof.
!  The larger lvlprof, the more accurate the interpolation
!  (we recommend using lvlprof=200 for a model run with about 50
!  points in the vertical direction).  Using the intermediate
!  sounding, the ARPS then generates a base state model sounding
!  for the particular vertical grid resolution chosen (i.e.,
!  the number of points in the vertical, nz, and the vertical grid
!  spacing, dz).
!
!  depthp:
!
!  The depth of atmosphere over which the interpolated profiles
!  will be defined.  depthp should be greater than or equal to
!  (nz-3)*dz, i.e., larger than the physical depth of the model
!  domain.  Otherwise, the code will extrapolate for gridpoints
!  outside the domain, leading to possible inconsistencies.  At all
!  costs, any such extrapolation should be avoided.
!
!-----------------------------------------------------------------------
!

  REAL :: uprof (lvlprof)      ! Temporary work array
  REAL :: vprof (lvlprof)      ! Temporary work array
  REAL :: ptprof(lvlprof)      ! Temporary work array
  REAL :: pprof (lvlprof)      ! Temporary work array
  REAL :: qvprof(lvlprof)      ! Temporary work array
  REAL :: zprof (lvlprof)      ! Temporary work array

  REAL :: temp1(lvlprof)       ! Temporary work array
  REAL :: temp2(lvlprof)       ! Temporary work array
  REAL :: temp3(lvlprof)       ! Temporary work array
  REAL :: temp4(lvlprof)       ! Temporary work array
  REAL :: temp5(lvlprof)       ! Temporary work array
  REAL :: temp6(lvlprof)       ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL    :: factor,deltaz,pt0,t0
  INTEGER :: i,j,k,kdata, istat, nunit, onvf, itrnmin, jtrnmin
  REAL    :: p0inv,cpdrd,ptbar0, nstatsq, ternmin, rho00, tem
  LOGICAL :: wk  ! Logical control for Weisman & Klemp sounding.
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
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  p0inv=1.0/p0
  cpdrd=1.0/rddcp

  IF( inibasopt == 5 ) THEN
    wk = .true.
  ELSE
    wk = .false.
  END IF

  onvf = 0
  CALL avgz(zp, onvf,                                                   &
            nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, zs)

  IF( inibasopt == 1 .OR. inibasopt == 5 ) THEN
!
!-----------------------------------------------------------------------
!
!  Initialize the base state from discrete sounding data.  Note that
!  The desired format of the sounding data is described inside
!  subroutine SOUNDG.
!
!-----------------------------------------------------------------------
!
    deltaz = depthp/(lvlprof-1)

    IF (myproc ==0) &
       WRITE(6,'(/'' THE INPUT SOUNDING DATA IS SAMPLED AT A '',           &
    &         f7.2,'' meter interval.''/)') deltaz
    
    DO k=1,lvlprof
      zprof(k) = (k-1)*deltaz
    END DO
!
!-----------------------------------------------------------------------
!
!  Obtain, from the input sounding, the base state vertical profiles
!  of the two horizontal wind components, potential temperature,
!  pressure, and water vapor specific humidity at vertical levels
!  defined by zprof.
!
!-----------------------------------------------------------------------
!
    CALL zprofil(zprof,lvlprof,                                         &
                 uprof,vprof,ptprof,pprof,qvprof,                       &
                 temp1,temp2,temp3,temp4,temp5,temp6,wk)
!
!
!-----------------------------------------------------------------------
!
!  Linearly interpolate the base state variables from the
!  intermediate evenly-spaced profile to the model grid.
!
!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!
!  Interpolation of scalars.
!
!-----------------------------------------------------------------------

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1

          kdata = INT( zs(i,j,k)/deltaz ) + 1
          kdata = MIN( lvlprof-1, MAX(1,kdata) )
          factor = (zs(i,j,k)-zprof(kdata))/deltaz

          ptbar(i,j,k)=ptprof(kdata)+                                   &
                       factor*(ptprof(kdata+1)-ptprof(kdata))

          pbar (i,j,k)=pprof(kdata)+                                    &
                       factor*(pprof(kdata+1)-pprof(kdata))

          qvbar(i,j,k)=qvprof(kdata)+                                   &
                       factor*(qvprof(kdata+1)-qvprof(kdata))

          pibar(i,j,k) = (pbar(i,j,k)*p0inv)**rddcp
          rhobar(i,j,k)=pbar(i,j,k)/(rd*ptbar(i,j,k)*pibar(i,j,k))
        END DO
      END DO
    END DO

    IF( inibasopt == 5 ) THEN 

!-----------------------------------------------------------------------
! Overwrite ptbar and qvbar by directly calculating them from 
! Weisman-Klemp analytical functions. Make sure that formulations and 
! parameters used as the same as in SOUNDG. Keep pbar since it's 
! more accurate being calculated from hydrostatic equation on high-resolution
! vertival grid used to set up the profiles inside soundg.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Create Weisman & Klemp (JAS 1982) thermodynamic profile over 20km. depth.
!-----------------------------------------------------------------------
!
      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            IF ( zs(i,j,k) <= htrop ) THEN
              ptbar(i,j,k) = ptground + (pttrop-ptground)*               &
                             ( (zs(i,j,k)/htrop)**(5.0/4.0) )
              qvbar(i,j,k) = 1.0 - 0.75*( (zs(i,j,k)/htrop)**(5.0/4.0) )
            ELSE
              ptbar(i,j,k) = pttrop*EXP( (zs(i,j,k)-htrop)*g/(ttrop*cp) )
              qvbar(i,j,k) = 0.25
            END IF

            pibar(i,j,k) = (pbar(i,j,k)*p0inv)**rddcp
            rhobar(i,j,k)=pbar(i,j,k)/(rd*ptbar(i,j,k)*pibar(i,j,k))
          END DO
        END DO
      END DO

      DO k = 2,nz-2 
        DO j = 1,ny-1
          DO i = 1,nx-1
            tem1(i,j,k) = ptbar(i,j,k) * pibar(i,j,k)
          END DO
        END DO
      END DO

      CALL getqvs(nx,ny,nz, 1,nx-1,1,ny-1,2,nz-2, pbar,tem1, tem2)

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            qvbar(i,j,k) = qvbar(i,j,k) * tem2(i,j,k)
            IF( zs(i,j,k) <= mixtop) THEN
              qvbar(i,j,k) = qvmixed
            ENDIF
            qvbar(i,j,k) = min( qvmixed, qvbar(i,j,k), tem2(i,j,k)*rhmixed )
          END DO
        END DO
      END DO
      

    ENDIF

!-----------------------------------------------------------------------
!
!  Set top and bottom boundary conditions
!
!  Zero gradient is imposed on ptbar and qvbar at the top and bottom
!  boundaries. pt and qv above the top and below the bottom boundary
!  are used only in the calculations of turbulent fluxes through
!  the boundary. Zero gradient conditions ensure zero heat and
!  misture fluxes through the boundary. When such fluxes need to be
!  included, they are handled by the surface physics parameterization.
!
!  Hydrostatic relation is used to obtain pbar and rhobar above the
!  top and below the bottom boundary.
!
!-----------------------------------------------------------------------

    DO j=1,ny-1
      DO i=1,nx-1
        ptbar(i,j,1)=ptbar(i,j,2)
        qvbar(i,j,1)=qvbar(i,j,2)
      END DO
    END DO

    DO i=1,nx-1
      DO j=1,ny-1
        pibar(i,j,1)=pibar(i,j,2)+g*(zs(i,j,2)-zs(i,j,1))               &
                    /(0.5*cp*(ptbar(i,j,2)+ptbar(i,j,1)) )
        pbar (i,j,1) = p0 * pibar(i,j,1)**cpdrd
        rhobar(i,j,1)= pbar(i,j,1)/(rd*ptbar(i,j,1)*pibar(i,j,1))
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        ptbar (i,j,nz-1)=ptbar(i,j,nz-2)
        qvbar (i,j,nz-1)=qvbar(i,j,nz-2)
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        pibar(i,j,nz-1)=pibar(i,j,nz-2)-g*(zs(i,j,nz-1)-zs(i,j,nz-2))   &
                       /(0.5*cp*(ptbar(i,j,nz-2)+ptbar(i,j,nz-1)) )
        pbar (i,j,nz-1) = p0 * pibar(i,j,nz-1)**cpdrd
        rhobar(i,j,nz-1)= pbar(i,j,nz-1)/                               &
                         (rd*ptbar(i,j,nz-1)*pibar(i,j,nz-1))
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  Interpolate u profile to model grid ubar
!  and subtract the domain translation speed when the option is on.
!
!-----------------------------------------------------------------------

    CALL avgsu(zs,nx,ny,nz, 1,ny-1, 1,nz-1, zuv, ubar) ! ubar used as temp

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx

          kdata = INT( zuv(i,j,k)/deltaz ) + 1
          kdata = MIN( lvlprof-1, MAX(1,kdata) )
          factor = (zuv(i,j,k)-zprof(kdata))/deltaz

          ubar(i,j,k)=uprof(kdata)+                                     &
                      factor*(uprof(kdata+1)-uprof(kdata))

        END DO
      END DO
    END DO

    IF(grdtrns == 1) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            ubar(i,j,k) = ubar(i,j,k) - umove
          END DO
        END DO
      END DO

    END IF

!-----------------------------------------------------------------------
!
!  Interpolate v profile to model grid ubar
!  and subtract the domain translation speed when the option is on.
!
!-----------------------------------------------------------------------

    CALL avgsv(zs,nx,ny,nz, 1,nx-1, 1,nz-1, zuv, vbar) ! vbar used as temp

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1

          kdata = INT( zuv(i,j,k)/deltaz ) + 1
          kdata = MIN( lvlprof-1, MAX(1,kdata) )
          factor = (zuv(i,j,k)-zprof(kdata))/deltaz

          vbar(i,j,k)=vprof(kdata)+                                     &
                      factor*(vprof(kdata+1)-vprof(kdata))

        END DO
      END DO
    END DO

    IF(grdtrns == 1) THEN

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            vbar(i,j,k) = vbar(i,j,k) - vmove
          END DO
        END DO
      END DO

    END IF


    DO j=1,ny-1
      DO i=1,nx
        ubar(i,j,1   )=ubar(i,j,2   )
        ubar(i,j,nz-1)=ubar(i,j,nz-2)
      END DO
    END DO

    DO j=1,ny
      DO i=1,nx-1
        vbar(i,j,   1)=vbar(i,j,   2)
        vbar(i,j,nz-1)=vbar(i,j,nz-2)
      END DO
    END DO

  ELSE IF( inibasopt == 2 ) THEN
!
!-----------------------------------------------------------------------
!
!  Specify a dry isentropic base state, using analytical formulations.
!  Hydrostatic relation is used.
!
!  This option can be switched on when needed.
!
!-----------------------------------------------------------------------
!
    pt0 = 300.0

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          pibar(i,j,k) = 1.0-g/(cp*pt0)*zs(i,j,k)
          pbar(i,j,k) = p0 * pibar(i,j,k)**cpdrd
          ptbar(i,j,k)=pt0
          rhobar(i,j,k)=pbar(i,j,k)/(rd*pt0*pibar(i,j,k))
          qvbar(i,j,k)=0.0
        END DO
      END DO
    END DO

  ELSE IF( inibasopt == 3) THEN
!
!-----------------------------------------------------------------------
!
!  Specify a dry isothermal base state, using analytical formulations.
!  Hydrostatic relation is used.
!
!  This option can be switched on when needed.
!
!-----------------------------------------------------------------------
!
    t0 = 250.0
    rho00 = p0/(rd*t0)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          pbar(i,j,k) = p0 * EXP(-g*zs(i,j,k)/(rd*t0))
          pibar(i,j,k)=(pbar(i,j,k)*p0inv)**rddcp
          ptbar(i,j,k)=t0/pibar(i,j,k)
          rhobar(i,j,k)=pbar(i,j,k)/(rd*t0)
          IF(bsnesq == 1) rhobar(i,j,k)=rho00
          qvbar(i,j,k)=0.0
        END DO
      END DO
    END DO

  ELSE IF(  inibasopt == 4 ) THEN ! Case of constant static stability

    nstatsq = 0.0001             ! constant static stability (1/s**2)
    ptbar0 = 300.0               ! PTBAR at z=0.0

    rho00 = p0/(rd*ptbar0)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          ptbar(i,j,k) = ptbar0*EXP(nstatsq*zs(i,j,k)/g )
        END DO
      END DO
    END DO

    tem = g*g/(nstatsq*cp*ptbar0)

    DO k=1,nz-1
      DO i=1,nx-1
        DO j=1,ny-1
          pibar(i,j,k)=1.0+tem*(EXP(-nstatsq/g*zs(i,j,k))-1.0)
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO i=1,nx-1
        DO j=1,ny-1
          pbar(i,j,k) =  p0*( pibar(i,j,k)**cpdrd )
          qvbar(i,j,k)=0.0
          rhobar(i,j,k)=pbar(i,j,k)/(rd*ptbar(i,j,k)*pibar(i,j,k))
          IF(bsnesq == 1) rhobar(i,j,k)=rho00
        END DO
      END DO
    END DO

  ELSE IF( inibasopt == 6 ) THEN

!-----------------------------------------------------------------------
!
!  Incompressible constant density and constant potential temperature
!  flow. Hydrostatic.
!
!-----------------------------------------------------------------------

    pt0 = 300.0

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          ptbar(i,j,k)=pt0
          rhobar(i,j,k)=p0/(rd*pt0*1.0)
          pbar(i,j,k)=p0 - rhobar(i,j,k)*g*zs(i,j,k)
          pibar(i,j,k) = (pbar(i,j,k)*p0inv)**rddcp
          qvbar(i,j,k)=0.0
        END DO
      END DO
    END DO

  ELSE

    IF (myproc ==0) WRITE(6,'(1x,a,i3,/1x,a)') &
        'Base state initialization option inibasopt=',inibasopt,        &
        'not implemented or available. Job stopped in INIBASE.'

    CALL arpsstop('arpsstop called from INIBASE with inibasopt',1)

  END IF

  IF ( inibasopt /= 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  For all cases except when external sounding (inibasopt=1) is used,
!  set the base state winds according to the viniopt option.
!
!-----------------------------------------------------------------------
!

    IF ( viniopt == 1 ) THEN
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            ubar(i,j,k) = ubar0
          END DO
        END DO
      END DO

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            vbar(i,j,k) = vbar0
          END DO
        END DO
      END DO

    ELSE IF ( viniopt == 2 ) THEN !  user specified wind profile
!
!-----------------------------------------------------------------------
!
!  When viniopt = 2, the profiles of ubar and vbar are supposed to be
!  specified by the user, by editing the following block of program.
!
!  As a default, ubar and vbar are set to zero.
!
!-----------------------------------------------------------------------
!
      CALL edgfill(zs,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
      CALL mpsendrecv2dew(zs,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(zs,nx,ny,nz,nbc,sbc,0,tem1)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            ubar(i,j,k) = 0.0
          END DO
        END DO
      END DO

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            vbar(i,j,k) = 0.0
          END DO
        END DO
      END DO

    ELSE IF ( viniopt == 3 ) THEN !  Weisman-Klemp wind profile
 
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=2,nx-1
            ubar(i,j,k) = ubar0*tanh((zs(i,j,k)+zs(i-1,j,k))*0.5/zshear)
          END DO
        END DO
      END DO

      DO k=1,nz-1
        DO j=2,ny-1
          DO i=1,nx-1
            vbar(i,j,k) = vbar0*tanh((zs(i,j,k)+zs(i,j-1,k))*0.5/zshear)
          END DO
        END DO
      END DO

      ubar(1, :,:)=ubar(2   ,:,:)
      ubar(nx,:,:)=ubar(nx-1,:,:)

      ubar(:,:,1   )=ubar(:,:,2   )
      ubar(:,:,nz-1)=ubar(:,:,nz-2)

      vbar(:,1, :)=vbar(:,2,   :)
      vbar(:,ny,:)=vbar(:,ny-1,:)

      vbar(:,:,   1)=vbar(:,:,   2)
      vbar(:,:,nz-1)=vbar(:,:,nz-2)

    ELSE

     IF (myproc ==0) WRITE (6, '(1x,a,i4,a)')  &
            'Unexpected option for viniopt: ',viniopt,                  &
            ', Model did not stopped in subroutine INIBASE.'

    END IF

    IF(grdtrns == 1) THEN

       ubar(:,:,:) = ubar(:,:,:) - umove
       vbar(:,:,:) = vbar(:,:,:) - vmove

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate rhostr etc.
!
!-----------------------------------------------------------------------
!
  DO k= 1,nz-1
    DO j= 1,ny-1
      DO i= 1,nx-1
        ptbari(i,j,k)=1.0/ptbar(i,j,k)
        pbari(i,j,k) =1.0/pbar(i,j,k)
        rhostr(i,j,k)=ABS(j3(i,j,k))*rhobar(i,j,k)
        rhostri(i,j,k)=1.0/rhostr(i,j,k)
      END DO
    END DO
  END DO

  CALL writesnd(nx,ny,nz,ubar,vbar,ptbar,pbar,qvbar,zp, rhobar, zs)

  RETURN
END SUBROUTINE inibase

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ZPROFIL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE zprofil(zprof,lvlprof,                                       &
           uprof,vprof,ptprof,pprof,qvprof,                             &
           zsnd,temp1,temp2,temp3,temp4,temp5,wk)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the initial sounding profile at regularly-spaced grid
!  levels.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    zprof     1-D array defining the regularly-spaced vertical grid
!              levels at which base state profiles will be sampled(m)
!    lvlprof   The number of grid levels in zprof.
!    WK        Logical flag for Weisman&Klemp analytical sounding
!
!  OUTPUT:
!
!    uprof     Profile of u-velocity defined at levels given
!              in zprof (m)
!    vprof     Profile of v-velocity defined at levels given
!              in zprof (m)
!    ptprof    Profile of potential temperature defined at levels
!              given in zprof (K)
!    pprof     Profile of pressure defined at levels given in zprof
!              (Pascal)
!    qvprof    Profile of water vapor specific humidity defined at
!              levels given in zprof (kg/kg)
!
!  WORK ARRAYS:
!
!    zsnd     Temporary work array to be used to store the actual
!             height of the input sounding data.
!    temp1    Temporary work array
!    temp2    Temporary work array
!    temp3    Temporary work array
!    temp4    Temporary work array
!    temp5    Temporary work array
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

  INTEGER :: lvlprof          ! The number of grid levels in zprof

  REAL :: zprof (lvlprof)     ! 1-D array defining the regularly-spaced
                              ! vertical grid levels at which base
                              ! state profiles will be sampled (m)

  REAL :: uprof (lvlprof)     ! Profile of u-velocity defined at
                              ! levels given in zprof (m)
  REAL :: vprof (lvlprof)     ! Profile of v-velocity defined at
                              ! levels given in zprof (m)
  REAL :: ptprof(lvlprof)     ! Profile of potential temperature
                              ! defined at levels given in zprof (K)
  REAL :: pprof(lvlprof)      ! Profile of pressure defined at
                              ! levels given in zprof (Pascal)
  REAL :: qvprof(lvlprof)     ! Profile of water vapor specific
                              ! humidity defined at levels given in
                              ! zprof (kg/kg)

  REAL :: zsnd (lvlprof)      ! Temporary work array
  REAL :: temp1(lvlprof)      ! Temporary work array
  REAL :: temp2(lvlprof)      ! Temporary work array
  REAL :: temp3(lvlprof)      ! Temporary work array
  REAL :: temp4(lvlprof)      ! Temporary work array
  REAL :: temp5(lvlprof)      ! Temporary work array

  LOGICAL :: wk               ! Flag for W&K analytical sounding
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: lvlsnd
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
!  Initialize the sounding profiles (analytical or from external
!  data).  All arguments are output.
!
!-----------------------------------------------------------------------
!
  CALL soundg(uprof,vprof,ptprof,pprof,qvprof,lvlprof,                  &
              zsnd,lvlsnd, temp1,temp2,temp3,temp4,temp5,wk )

!
!-----------------------------------------------------------------------
!
!  Interpolate the sounding profile at levels zsnd to zprof.
!  Outputs are also stored in arrays *prof, where * = u,v,w, etc.
!
!-----------------------------------------------------------------------
!
  CALL sndintrp(uprof,vprof,ptprof,pprof,qvprof,zsnd,lvlsnd,            &
                zprof,lvlprof, temp1)

  RETURN
END SUBROUTINE zprofil

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOUNDG                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE soundg(usnd,vsnd,ptsnd,psnd,qvsnd,lvlprof,                   &
           zsnd,lvlsnd, tsnd,rhsnd,ptvsnd,pisnd,dpsnd,wk)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Specify the sounding profiles for u, v, pt, p and qv from analytical
!  functions or from an external sounding file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  10/31/94 (E. Adlerman)
!  Added option for dewpoint in input sounding
!  Added option for Weisman & Klemp sounding
!
!  9/22/1995 (M. xue)
!  Clarified the definition of sounding data recard 7.
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    lvlprof  The dimension of arrays usnd,vsnd,ptsnd,psnd,qvsnd,zsnd
!             This parameter is the number of levels in the
!             intermediate interpolated sounding.
!    WK       Logical flag for Weisman&Klemp thermodynamical
!             profile. Winds are set to zero.
!  OUTPUT:
!
!    usnd     u-velocity sounding data at original data levels (m/s)
!    vsnd     v-velocity sounding data at original data levels (m/s)
!    ptsnd    Potential temperature sounding data at original data
!             levels (K)
!    psnd     Pressure sounding data at original data levels (Pascal)
!    qvsnd    Water vapor specific humidity sounding data at the
!             original data levels (kg/kg)
!    zsnd     The height (above ground) of the sounding data levels
!             (m)
!    lvlsnd   The number of vertical levels in the input sounding
!
!  WORK ARRAYS:
!
!    tsnd     Temporary work array
!    ptvsnd   Temporary work array
!    pisnd    Temporary work array
!    rhsnd    Temporary work array
!    tdsnd    Temporary work array
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

  INCLUDE 'mp.inc'          ! Message passing parameters.

  INTEGER :: lvlprof

  REAL :: usnd (lvlprof)     ! u-velocity from sounding data (m/s)
  REAL :: vsnd (lvlprof)     ! v-velocity from sounding data (m/s)
  REAL :: ptsnd(lvlprof)     ! Potential temperature sounding data (K)
  REAL :: psnd (lvlprof)     ! Pressure from sounding data (Pascal)
  REAL :: qvsnd(lvlprof)     ! Water vapor specific humidity
                             ! from sounding data (kg/kg)
  REAL :: zsnd (lvlprof)     ! The height of the sounding data levels
                             ! (m)

  REAL :: tsnd  (lvlprof)    ! Temporary work array
  REAL :: ptvsnd(lvlprof)    ! Temporary work array
  REAL :: pisnd (lvlprof)    ! Temporary work array
  REAL :: rhsnd (lvlprof)    ! Temporary work array
  REAL :: dpsnd (lvlprof)    ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
!
!  The following three strings designate the type of sounding data.
!
!  if height(1:1)='h' or 'H', the sounding is given on height levels
!  if height(1:1)='p' or 'P', the sounding is given on pressure levels
!
!  if therm(1:1) ='t' or 'T', the sounding is specified in temperature
!  if therm(1:1) ='p' or 'P', the sounding is specified in potential
!                             temperature.
!
!  if humid(1:1) ='s' or 'S', the soundings uses specific humidity
!  if humid(1:1) ='r' or 'R', the sounding uses relative humidity,
!  if humid(1:1) ='d' or 'D', the soundings uses dewpoint temperature,
!
!  if wind(1:1) ='x' or 'Y', the sounding is specified in x and y
!                            component of velocity (u and v)
!  if wind(1:1) ='d' or 'D', the sounding is specified in direction
!                            and speed in m/s.
!  if wind(1:1) ='k' or 'K', the sounding is specified in direction
!                            and speed in knots.
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=30) :: height, therm, humid , wind
  LOGICAL :: in_is_p, in_is_t, in_is_rh, in_is_dp, wk
  CHARACTER (LEN=80) :: dummy
  CHARACTER (LEN=72) :: stime, sdate, sloc
  INTEGER :: k,iter,niter,lvlsnd,istat
  REAL    :: psnd1,zsnd1
  INTEGER :: lenstr
  REAL    :: p0inv,cpdrd
  INTEGER :: item1, item2
  REAL    :: tem1
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
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  p0inv=1.0/p0
  cpdrd=1.0/rddcp
!
!-----------------------------------------------------------------------
!
!  Create Weisman & Klemp thermodynamic profile over 20km. depth.
!  Winds are set to zero. See W&K, JAS 1982 for details.
!
!-----------------------------------------------------------------------
!

  IF ( wk ) THEN

    lvlsnd = 201
    zsnd1  = 0.0
    psnd1  = 100000.0
    height = 'height'
    humid  = 'rh'
    therm  = 'pt'
    wind   = 'uv'

!   pttrop = 343.0           ! Tropopause pot.temp.
!   ttrop  = 213.0           ! Tropopause temp.
!   ptground = 300.0         ! Surface pot.temp.
!   htrop  = 12000.0         ! Tropopause height
!   qvmixed= 0.015           ! Mixed layer mixing ratio
!   mixtop = 1200.0          ! Mixed layer height


    zsnd(1) = zsnd1
    DO k=2, lvlsnd
      zsnd(k)=(20000.0/REAL(lvlsnd-1)) * REAL(k-1)
    END DO

    DO k=1, lvlsnd
      IF ( zsnd(k) <= htrop ) THEN
        ptsnd(k) =                                                      &
            ptground + (pttrop-ptground)*( (zsnd(k)/htrop)**(5.0/4.0) )
        qvsnd(k) =                                                      &
            1.0 - 0.75*( (zsnd(k)/htrop)**(5.0/4.0) )
      ELSE
        ptsnd(k) = pttrop*EXP( (zsnd(k)-htrop)*g/(ttrop*cp) )
        qvsnd(k) = 0.25
      END IF
    END DO

    DO k=1, lvlsnd
      usnd(k) = ubar0*tanh(zsnd(k)/zshear)  ! not really used
      vsnd(k) = vbar0*tanh(zsnd(k)/zshear)  ! not really used
    END DO

  ELSE

!
!-----------------------------------------------------------------------
!
!              Sounding File Format for ARPS 3.0
!
!  Record 1: a header line, such as "1-D Sounding Input for ARPS"
!            (skipped)
!  Record 2: miscellaneous description of sounding (skipped)
!  Record 3: time of sounding (character*72)
!  Record 4: date of sounding (character*72)
!  Record 5: location of sounding (character*72)
!  Record 6: three character strings designate the sounding
!            data type, e.g.
!            'pressure' 'potential_temperature' 'relative_humidity'.
!            Only the first character of the strings
!            is decoded, and thus the first character should not be
!            left blank.  Note that either upper or lower case may be
!            used.  A more detailed explanation is provided in the
!            portion of the code where these strings are declared.
!  Record 7: Ground-level height (m) and the correspoding pressure
!            (Pascal) when sounding is specified at height levels.
!            When it is given on pressure levels, this record is
!            not used. The last sounding data is assumed to be at
!            the ground level (z=0).
!  Record 8: number of data levels in sounding (variable lvlsnd )
!  Record 9: a line of data column labels (not read in)
!  Record 10 to the end:
!            sounding data in the order 'z/p, pt/t, qv/rh, u, v'
!
!  For records 10  to the end, there is one line of data
!            corresponding to each sounding level.
!
!  Important convention:
!            Line 10 corresponds to the level k = lvlsnd
!            Line 11 corresponds to the level k = (lvlsnd -1)
!            etc.
!
!  Units of the data:
!            pressure: Pascals, height: meters, temperature, dewpoint,
!            and potential temperature: degrees Kelvin, specific
!            humidity kg/kg, relative humidity: 0 to 1.
!
!  CAUTION: The sounding data have to be listed in the order of
!           decreasing height or increasing pressure.
!
!-----------------------------------------------------------------------
!

    istat = 0

!
!-----------------------------------------------------------------------
!
!  Open the sounding file
!
!-----------------------------------------------------------------------
!
    CALL getunit( sfunit )

    OPEN(UNIT=sfunit, FILE=trim(sndfile), STATUS='old', IOSTAT=istat)

    IF(istat /= 0) THEN
      IF (myproc == 0 ) THEN
        WRITE (6,1000) 'Error in opening sounding file'
        WRITE (6,'(1x,a,a)') 'Check file ',sndfile
      END IF
      CALL arpsstop('arpsstop called from INIBASE with opening snd file',1)
    ELSE
      lenstr = 80
      CALL strlnth( sndfile, lenstr )
      IF (myproc == 0 ) &
         WRITE(6,'(1x,/2a,a,i3/)') 'Sounding file ', sndfile(1:lenstr),    &
              ' opened using FORTRAN unit',sfunit
    END IF

    READ (sfunit,1000,IOSTAT=istat) dummy ! header line
    IF (myproc ==0) WRITE(6,1000)dummy

    READ (sfunit,1000,IOSTAT=istat) dummy ! miscellaneous descriptions
    IF (myproc ==0) WRITE(6,1000)dummy

    READ (sfunit,1000,IOSTAT=istat) stime ! time of sounding
    IF (myproc ==0) WRITE(6,1000)stime


    READ (sfunit,1000,IOSTAT=istat) sdate ! date of sounding
    IF (myproc ==0) WRITE(6,1000)sdate

    READ (sfunit,1000,IOSTAT=istat) sloc  ! location of sounding
    IF (myproc ==0) WRITE(6,1000)sloc
!
!-----------------------------------------------------------------------
!
!  The following is for maintaining backwards compatibility with
!  ARPS 3.* sounding files which did not have a wind data type
!  indicator.
!
!-----------------------------------------------------------------------
!
    READ (sfunit,'(a)',IOSTAT=istat) dummy   ! tem string for data types

    item1 = 0
    item1 = INDEX( dummy(item1+1:80), '''' ) + item1
    item2 = INDEX( dummy(item1+1:80), '''' ) + item1
    height = dummy((item1+1):item2-1)

    item1 = item2
    item1 = INDEX( dummy(item1+1:80), '''' ) + item1
    item2 = INDEX( dummy(item1+1:80), '''' ) + item1
    therm = dummy(item1+1:item2-1)

    item1 = item2
    item1 = INDEX( dummy(item1+1:80), '''' ) + item1
    item2 = INDEX( dummy(item1+1:80), '''' ) + item1
    humid = dummy((item1+1):item2-1)

    item1 = item2
    item1 = INDEX( dummy(item1+1:80), '''' ) + item1
    item2 = INDEX( dummy(item1+1:80), '''' ) + item1
    IF ( item1 == item2 ) THEN
      wind = 'uv'
    ELSE
      wind  = dummy(item1+1:item2-1)
    END IF

    IF (myproc == 0 ) WRITE(6,'(1x,4a)') height, therm, humid, wind

    READ (sfunit,*,IOSTAT=istat) zsnd1,psnd1
    IF (myproc ==0)  &
       WRITE(6,'(1x,''SNDL='',F13.2,''PSNDL='',F13.2)') zsnd1,psnd1

    READ (sfunit,*,IOSTAT=istat) lvlsnd   ! actual number of data
                                          ! levels in the sounding
                                          ! file
    WRITE(6,'(1x,i6)')

    IF(lvlsnd > lvlprof) THEN

      IF (myproc == 0 ) THEN
        WRITE (6,1000) 'Error in sounding file'
        WRITE (6,'(1x,a,i6)') 'lvlsnd cannot be greater than ',lvlprof
        WRITE (6,1000) 'Check subroutine SOUNDG'
      END IF

      CALL arpsstop('arpsstop called from INIBASE with number of levels &
        & in the sounding ',1)

    END IF


    READ (sfunit, 1000, IOSTAT=istat) dummy ! a separator line,
    IF (myproc ==0) WRITE(6,1000)dummy

    IF(istat /= 0) THEN

      IF (myproc == 0 ) THEN
        WRITE (6,1000) 'Error in sounding file'
        WRITE (6,1000) 'Error in sounding file header'
        WRITE (6,1000) 'Check sounding file', sndfile
      END IF

      CALL arpsstop('arpsstop called from INIBASE with opening snd file',1)

    END IF

!-----------------------------------------------------------------------
!
!  Read in the sounding data:
!
!-----------------------------------------------------------------------

    DO k = lvlsnd ,1,-1

      READ (sfunit, *, IOSTAT=istat)                                    &
            zsnd(k), ptsnd(k), qvsnd(k), usnd(k), vsnd(k)

      IF (myproc == 0) WRITE(6,'(1x,i4,5f13.6)')  &
          k,zsnd(k), ptsnd(k), qvsnd(k), usnd(k), vsnd(k)

      IF(istat /= 0) THEN

        IF (myproc == 0) THEN
          WRITE (6,1000) 'Error in sounding file'
          WRITE (6,1000) 'Error in sounding file data'
          WRITE (6,1000) 'Check ',sndfile
        END IF
        CALL arpsstop('arpsstop called from INIBASE opening snd file',1)

      END IF

    END DO

    CLOSE (UNIT=sfunit)      ! close the sounding file
    CALL retunit( sfunit )

    IF (myproc == 0) WRITE(6,'(/1x,a)')'Sounding data reading complete.'
!
!-----------------------------------------------------------------------
!
!  End of sounding data input.
!
!-----------------------------------------------------------------------
!

  END IF
!
!-----------------------------------------------------------------------
!
!  If wind input is direction and speed, then find u,v.
!  For now there is no rotation of winds due to the sounding
!  not coinciding with the map projection longitude.
!  We tell the conversion program the sounding is at trulat1
!  and trulon.  In the future this can be the actual sounding
!  lat, lon or the center of the ARPS grid.
!  Here tsnd, pisnd and ptvsnd are used as dummy work arrays.
!
!-----------------------------------------------------------------------
!
  IF( wind(1:1) == 'd' .OR. wind(1:1) == 'D' .OR.                       &
        wind(1:1) == 'k' .OR. wind(1:1) == 'K' ) THEN

    DO k=1,lvlsnd
      pisnd(k)=trulon
    END DO

    IF( wind(1:1) == 'k' .OR. wind(1:1) == 'K' ) THEN
      DO k=1,lvlsnd
        vsnd(k) = vsnd(k)*0.514444
      END DO
    END IF

    CALL ddrotuv(lvlsnd,pisnd,usnd,vsnd,ptvsnd,                         &
                 usnd,vsnd)

  END IF
!
!-----------------------------------------------------------------------
!
!  If input data are not in height units, copy into the pressure
!  array and set the input data type flag.
!
!-----------------------------------------------------------------------
!

  IF( height(1:1) /= 'h' .AND. height(1:1) /= 'H' ) THEN

    in_is_p = .true.
    DO k=1,lvlsnd
      psnd(k) =zsnd(k)
    END DO

  ELSE

    in_is_p = .false.

  END IF
!
!-----------------------------------------------------------------------
!
!  If input is temperature, copy it from the potential temperature
!  array and set the input data type flag.
!
!-----------------------------------------------------------------------
!
  IF( therm(1:1) == 't' .OR. therm(1:1) == 'T' ) THEN

    in_is_t = .true.
    DO k=1,lvlsnd
      tsnd(k) =ptsnd(k)
    END DO

  ELSE

    in_is_t = .false.

  END IF

!
!-----------------------------------------------------------------------
!
!  If input is relative humidity, copy it from the specific humidity
!  array and set the input data type flag. If input is dewpoint, copy
!  it from the specific humidity array and set the input data type
!  flag. Otherwise, humidity is already mixing ratio/specific humid.
!  and the data flag types are set appropriately.
!
!-----------------------------------------------------------------------
!
  IF( humid(1:1) == 'r' .OR. humid(1:1) == 'R' ) THEN


    in_is_rh = .true.
    in_is_dp = .false.
    DO k=1,lvlsnd
      rhsnd(k) =qvsnd(k)
    END DO

  ELSE IF ( humid(1:1) == 'd' .OR. humid(1:1) == 'D' ) THEN

    in_is_dp = .true.
    in_is_rh = .false.
    DO k=1, lvlsnd
      dpsnd(k) = qvsnd(k)
    END DO

  ELSE

    in_is_dp = .false.
    in_is_rh = .false.
    DO k=1,lvlsnd
      rhsnd(k) = 0.0
      dpsnd(k) = 0.0
    END DO


  END IF

!
!-----------------------------------------------------------------------
!
!  Calculate the derived variables from the input sounding
!
!-----------------------------------------------------------------------
!

  IF( in_is_p ) THEN
!
!-----------------------------------------------------------------------
!
!  If the sounding data are specified on pressure levels...
!
!-----------------------------------------------------------------------
!
    DO k=1,lvlsnd

      IF( psnd(k) <= 0.0) THEN

        IF (myproc == 0) WRITE(6,'(2(/1x,a))')  &
             'Negative pressure found in the sounding profile.',        &
             'Job stopped. Check subroutine INIBASE.'
        CALL arpsstop('arpsstop called from INIBASE due to neg. pressure ',1)

      END IF

    END DO

    IF( in_is_t ) THEN

      DO k=1,lvlsnd
        ptsnd(k) = tsnd(k)*(p0/psnd(k))**rddcp
      END DO

    ELSE

      DO k=1,lvlsnd
        tsnd(k) = ptsnd(k)*(psnd(k)*p0inv)**rddcp
      END DO

    END IF

    IF( in_is_rh ) THEN

      CALL getqvs(1,1,lvlsnd, 1,1,1,1,1,lvlsnd, psnd,tsnd,qvsnd)

      DO k=1,lvlsnd
        qvsnd(k) = qvsnd(k) * rhsnd(k)
      END DO

    ELSE IF ( in_is_dp ) THEN

      CALL getqvs(1,1,lvlsnd, 1,1,1,1,1,lvlsnd, psnd,dpsnd,qvsnd)

    END IF
!
!-----------------------------------------------------------------------
!
!  The virtual temperature Tv:
!
!-----------------------------------------------------------------------
!
    DO k = 1,lvlsnd
      ptvsnd(k) = tsnd(k)*(1.0+rvdrd*qvsnd(k))/(1.0+qvsnd(k))
    END DO

!
!-----------------------------------------------------------------------
!
!  Integrate hydrostatic relation dz/dp = - RTv/(pg) for z
!
!-----------------------------------------------------------------------
!
    zsnd(1) = zsnd1

    DO k = 2,lvlsnd
      zsnd(k)=zsnd(k-1) - rd/g                                          &
            *(ptvsnd(k)+ptvsnd(k-1))/(psnd(k)+psnd(k-1))                &
            *(psnd(k)-psnd(k-1))
    END DO
!

  ELSE IF( .NOT. in_is_p ) THEN

!
!-----------------------------------------------------------------------
!
!  If the sounding data are specified at height levels
!
!-----------------------------------------------------------------------
!
    IF( in_is_rh .OR. in_is_dp ) THEN
!
!-----------------------------------------------------------------------
!
!  If the input fields are height and relative humidity/dewpoint, the
!  pressure, temperature and specific humidity profiles are coupled;
!  thus, one must iterate to compute the pressure and temperature.
!  We use a predetermined number (10) of iterations, which should
!  ensure convergence in even the most difficult situations.
!
!-----------------------------------------------------------------------
!
      niter = 10

      DO k=1,lvlsnd
        qvsnd(k) =0.0
      END DO
    ELSE
      niter = 1
    END IF


    DO iter = 1, niter

      IF( .NOT. in_is_t) THEN
!
!-----------------------------------------------------------------------
!
!  Thermodynamic variable is potential temperature...
!
!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!
!  Calculate ptvsnd, the virtual potential temperature sounding array
!
!-----------------------------------------------------------------------
!
        DO k = 1,lvlsnd
          ptvsnd(k) = ptsnd(k)*(1.0+rvdrd*qvsnd(k))/(1.0+qvsnd(k))
        END DO

!
!-----------------------------------------------------------------------
!
!  Integrate the hydrostatic relation to get PISND, the
!  nondimensional pressure: d(pi)/dz = -g/(cp * ptv)
!
!-----------------------------------------------------------------------
!
!

        pisnd(1) = (psnd1*p0inv)**rddcp

        DO k = 2,lvlsnd
          pisnd(k) = pisnd(k-1) - g/cp                                  &
                     *(zsnd(k)-zsnd(k-1))                               &
                     /(0.5*(ptvsnd(k)+ptvsnd(k-1)))
        END DO

        DO k = 1,lvlsnd
          psnd(k) = p0 * (pisnd(k)**cpdrd)
          tsnd(k) = ptsnd(k) * pisnd(k)
        END DO

      ELSE
!
!-----------------------------------------------------------------------
!
!  Thermodynamic variable is air temperature...
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Calculate the virtual temperature and store it in ptvsnd
!
!-----------------------------------------------------------------------
!
        DO k = 1,lvlsnd
          ptvsnd(k) = tsnd(k)*(1.0+rvdrd*qvsnd(k))/(1.0+qvsnd(k))
        END DO

!
!-----------------------------------------------------------------------
!
!  Calculate log(p) according to  d(log(p))/dz = -g/(r * tv)
!
!-----------------------------------------------------------------------
!
        pisnd(1) = LOG( psnd1 )

        DO k = 2,lvlsnd
          pisnd(k) = pisnd(k-1) - g/rd                                  &
                     *(zsnd(k)-zsnd(k-1))                               &
                     /(0.5*(ptvsnd(k)+ptvsnd(k-1)))
        END DO

        DO k = 1,lvlsnd
          psnd(k) = EXP( pisnd(k) )
          ptsnd(k) = tsnd(k) * (p0/psnd(k))**rddcp
        END DO

      END IF

      IF ( in_is_rh ) THEN

        CALL getqvs(1,1,lvlsnd, 1,1,1,1,1,lvlsnd, psnd,tsnd,qvsnd)

        DO k=1,lvlsnd
          qvsnd(k) = qvsnd(k) * rhsnd(k)
        END DO

      ELSE IF ( in_is_dp ) THEN

        CALL getqvs(1,1,lvlsnd, 1,1,1,1,1,lvlsnd, psnd,dpsnd,qvsnd)

      END IF

    END DO

  END IF

  IF( .NOT. in_is_rh ) THEN

    CALL getqvs(1,1,lvlsnd, 1,1,1,1,1,lvlsnd, psnd,tsnd,rhsnd)

    DO k=1,lvlsnd
      rhsnd(k) = qvsnd(k)/rhsnd(k)
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Add fixed mixed layer to WK profile and recalculate RH values
!
!-----------------------------------------------------------------------
!
  IF ( wk ) THEN

    DO k=1, lvlsnd

      IF (zsnd(k) <= mixtop) THEN
        tem1     = f_qvsat(psnd(k),tsnd(k))  ! saturated specific humidity
        qvsnd(k) = min(qvmixed, tem1 * rhmixed)
        rhsnd(k) = qvsnd(k)/ tem1
      END IF

    END DO

  END IF


!
!-----------------------------------------------------------------------
!
!  Write out the uninterpolated sounding profile
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    WRITE(6,'(/1x,a/)')                                                 &
      'Base state profiles at original sounding data levels:'

    WRITE(6,'(1x,2a)')                                                  &
       '   n    z(m)   p(Pascal)  PT(K)   T(K)',                        &
       '    Qv(kg/kg)   RH   ubar(m/s) vbar(m/s)'

    DO k=lvlsnd,1,-1
      WRITE(6,'(1x,i4,f9.2,f9.1,2f9.3,f9.6,3f9.3)')                     &
          k,zsnd(k),psnd(k),ptsnd(k),                                   &
          tsnd(k), qvsnd(k), rhsnd(k), usnd(k), vsnd(k)
    END DO
  END IF

  RETURN

  1000  FORMAT(1X,a)

END SUBROUTINE soundg

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SNDINTRP                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sndintrp(uprof,vprof,ptprof,pprof,qvprof,zsnd,lvlsnd,        &
           zprof,lvlprof, temp1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate the input sounding profile on the grid "zsnd" to a
!  profile on the model grid "zprof".
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    uprof    u-velocity sounding data at original data levels (m/s)
!    vprof    v-velocity sounding data at original data levels (m/s)
!    ptprof   Potential temperature sounding data at original data
!             levels (K)
!    pprof    Pressure sounding data at original data levels (Pascal)
!    qvprof   Water vapor specific humidity sounding data at the
!             original data levels (kg/kg)
!    zsnd     The height of the input sounding data levels (m)
!    lvlsnd   The number of input sounding data levels
!
!    zprof    The height of equally spaced grid levels to which
!             the sounding data are interpolated (m).
!    lvlprof  The number of such grid levels.
!
!  OUTPUT:
!
!    uprof    u-velocity data defined on grid levels given in zprof
!    vprof    v-velocity data defined on grid levels given in zprof
!    ptprof   Potential temperature data defined on grid levels
!             given in zprof
!    pprof    Pressure data defined on grid levels given in zprof
!    qvprof   Water vapor specific humidity data defined on grid
!             levels given in zprof
!
!  WORK ARRAY:
!
!    temp1     Temporary work array
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

  INTEGER :: lvlprof, lvlsnd

  REAL :: uprof (lvlprof)     ! Profile of u-velocity defined at
                              ! levels given in zprof (m)
  REAL :: vprof (lvlprof)     ! Profile of v-velocity defined at
                              ! levels given in zprof (m)
  REAL :: ptprof(lvlprof)     ! Profile of potential temperature
                              ! defined at levels given in zprof (K)
  REAL :: pprof(lvlprof)      ! Profile of pressure defined at
                              ! levels given in zprof (Pascal)
  REAL :: qvprof(lvlprof)     ! Profile of water vapor specific
                              ! humidity defined at levels given in
                              ! zprof (kg/kg)
  REAL :: zsnd (lvlprof)      ! Temporary work array
  REAL :: zprof (lvlprof)     ! 1-D array defining the regularly-spaced
                              ! vertical grid levels at which base
                              ! state
                              ! profiles will be sampled (m)
  REAL :: temp1 (lvlsnd)      ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: k
!
!-----------------------------------------------------------------------
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,lvlsnd
    temp1(k) = uprof(k)
  END DO

  CALL inte1d(temp1,zsnd,lvlsnd, uprof,zprof,lvlprof)

  DO k=1,lvlsnd
    temp1(k) = vprof(k)
  END DO

  CALL inte1d(temp1,zsnd,lvlsnd, vprof,zprof,lvlprof)

  DO k=1,lvlsnd
    temp1(k) = ptprof(k)
  END DO

  CALL inte1d(temp1,zsnd,lvlsnd, ptprof,zprof,lvlprof)

  DO k=1,lvlsnd
    temp1(k) = pprof(k)
  END DO

  CALL inte1d(temp1,zsnd,lvlsnd, pprof,zprof,lvlprof)

  DO k=1,lvlsnd
    temp1(k) = qvprof(k)
  END DO

  CALL inte1d(temp1,zsnd,lvlsnd, qvprof,zprof,lvlprof)

  RETURN
END SUBROUTINE sndintrp
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INTE1D                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE inte1d(fdat,zdat,ndat,fpro,zpro,npro)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate data from sounding to a model grid-column with uniform
!  grid spacing.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  Input :
!
!    fdat     1-D array defined at levels in zdat to be
!             interpolated to levels defined by zpro.
!    zdat     The height of the input data given by fdat.
!    ndat     The number of data in array fdat.
!
!    zpro     The grid level height to which data are interpolated
!    npro     The number of interpolated data levels
!
!  Output:
!
!    fpro     Output data interpolated to grid levels defined by zpro.
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

  INTEGER :: ndat             ! The number of data in array fdat
  INTEGER :: npro             ! The number of interpolated data levels

  REAL :: fdat(ndat)          ! 1-D array defined at levels in zdat to
                              ! be interpolated to levels defined by
                              ! zpro
  REAL :: zdat(ndat)          ! The height of the input data given by
                              ! fdat
  REAL :: zpro(npro)          ! The grid level height to which data are
                              ! interpolated
  REAL :: fpro(npro)          ! Output data interpolated to grid levels
                              ! defined by zpro
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i0,i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  i0=1
  DO j=1,npro

    IF(zpro(j) <= zdat(1))THEN
      fpro(j)=fdat(1)
      CYCLE
    END IF

    IF(zpro(j) >= zdat(ndat))THEN
      fpro(j)=fdat(ndat)
      CYCLE
    END IF

!
!-----------------------------------------------------------------------
!
!  Find index for interpolation
!
!-----------------------------------------------------------------------
!
    DO i=i0,ndat-1
      IF(zpro(j) >= zdat(i).AND.zpro(j) < zdat(i+1))GO TO 15
    END DO

    CYCLE
    15      i0=i

    fpro(j)=fdat(i)+(fdat(i+1)-fdat(i))*(zpro(j)-zdat(i))               &
            /(zdat(i+1)-zdat(i))

  END DO

  RETURN
END SUBROUTINE inte1d
