!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SFCPHYSICS                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sfcphysics(nx,ny,nz,nzsoil,nstyps,                           &
           u,v,w,ptprt,pprt,qv,                                         &
           rhostr,ptbar,pbar,qvbar,                                     &
           x,y,z, zp,zpsoil,j1,j2,j3,j3soil,j3inv,j3soilinv,prcrate,    &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           radsw,rnflx,radswnet,radlwin,                                &
           cdm,cdh,cdq,usflx,vsflx,ptsflx,qvsflx,pbldpth,tsdiffus,      &
           tem1soil,tem2soil,tem3soil,tem4soil,tem1,tem2,tem3,tem4,     &
           tem5,shflx,lhflx,gflx,ct,evaprg,evaprtr,evaprr,qvsat,qvsata, &
           f34,temxy1,temxy2,temxy3,temxy4,temxpy)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the surface momentum, heat and surface moisture fluxes.
!  These fluxes will be passed into turbulent mixing subroutines.
!  The drag coefficients are returned in cdm, cdh and cdq.
!
!  The soil-vegetation model is also called to update its state
!  variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/07/1992.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  6/4/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  9/14/1992 (M. Xue)
!  Different drag coefficients defined for momentum, temperature and
!  moisture. Definition of ptsfc0 and qvsfc0 changed. qvbar added to
!  the argument list.
!
!  9/28/1992 (M. Xue)
!  Sign correction in all the flux formulations.
!
!  8/28/1993 (M. Xue and Hao jin)
!  Modified the flux formulations considering the terrain effect.
!
!  01/24/94 (Y. Liu and V. Wong)
!  Added the surface energy and moisture budget model to predict the
!  ground temperature and moistures.
!
!  9/04/1994 (K. Brewster, V. Wong and X. Song)
!  Facelift.
!
!  11/01/1994 (Y. Liu)
!  Subroutine name of soil model were changed from SFCEBM to SOILEBM
!  and the arguments passed were changed to 2-D arrays instead of 3-D
!  temporary arrays.
!
!  12/07/1994 (Y. Liu)
!  Cleaned up internal documentation and an empty DO loop of labeled
!  360 in subroutine SFCFLXSD.
!
!  01/28/1995 (V. Wong and X. Song)
!  Add the option of stability and roughness dependent surface layer
!  parameterization for both land and sea surfaces.
!
!  02/07/1995 (Y. Liu)
!  Fixed a bug in the calculation of surface precipitation,
!  tem1(i,j,4).
!
!  02/27/95 (V. Wong, Y. Liu and X. Song)
!  Delete subroutines of calculating C_u and C_pt at the neutral state
!  Make predicting procedure more consistent. Add a limiting parameter
!  "blimit" to prevent the drag coefficients from becoming too small
!  in the stable region.
!
!  03/02/1995 (Y. Liu)
!  Added the minimum surface total wind speed to guarantee the basic
!  heat and moisture fluxes.
!
!  02/07/96 (V.Wong and X.Song)
!  Changed the calculation of roughness zo over the sea according to
!  Clarke(1970).  The new formulation is insensitive to the value of
!  the depth. of the first model layer above the ground.
!  Added kvwtr to denote the Von Karman constant over the sea.
!  Set a lower limiter, zolimit, for zo, and an upper limiter, z1limit,
!  for depth of the surface layer z1.
!
!  12/08/98 (Donghai Wang and Vince Wang)
!  Added the snow cover.
!
!  2/25/1999 (M. Xue)
!  Reorganized this subroutine from original SFCFLX.
!  Calculations of drag coefficients and fluxes are grouped into
!  into new subroutines DRAGCOEF and SFCFLX, respectively.
!  Named local arrays are defined for many 2D variables for better
!  readability.
!
!  2000/02/29 (Gene Bassett, Yang Kun, Vince Wong)
!  Corrected typo: "c8=gammahl*stabp" from equation (15) in Byun (1990).
!
!  2001/12/05 (Yunheng Wang)
!  Merged M. Xue's version to the official release.
!
!  2001/12/10 (Ming Xue)
!  Streamlined the usage of work arrays.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        z component of velocity at a given time level (m/s)
!
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    pbldpth  Planetary boundary layer depth (m) at a given time level
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        The x-coord. of the physical and computational grid.
!             Defined at u-point.
!    y        The y-coord. of the physical and computational grid.
!             Defined at v-point.
!    zp       The physical height coordinate defined at w-point
!    zpsoil   The physical height coordinate defined at w-point
!
!    j1       Coordinate transformation Jacobian -d(zp)/d(x)
!    j2       Coordinate transformation Jacobian -d(zp)/d(y)
!    j3       Coordinate transformation Jacobian, d(zp)/d(z)
!    j3soil   Coordinate transformation Jacobian, d(zpsoil)/d(z)
!    j3inv        Inverse of j3
!    j3soilinv    Inverse of j3soil
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by the surface
!    radswnet Net solar radiation, SWin - SWout
!    radlwin  Incoming longwave radiation
!
!  OUTPUT:
!
!    cdm      Drag coefficient for momentum defined as scalar point
!    cdh      Drag coefficient for heat
!    cdq      Drag coefficient for moisture
!    usflx    Surface flux of u-momentum
!    vsflx    Surface flux of v-momentum
!    ptsflx   Surface flux of heat (K*kg/(m**2*s))
!    qvsflx   Surface flux of moisture (K*kg/(m**2*s))
!    pbldpth  Planetary boundary layer depth (m)
!
!  INPUT & OUTPUT:
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy moisture
!    snowdpth Snow depth (m)
!    qvsfc    Effective specific humidity at sfc.
!
!  WORK ARRAY:
!
!    ptsfc    Potential temperature at the ground level (K)
!    psfc     Surface pressure (Pascal)
!    wdspsfc  Wind speed at the surface
!    rhosfc   Base-state air density at the surface
!    pt1      Potential temperature at the first model level AGL
!    tair     Temperature of surface air
!    qvair    Specific humidity of surface air
!
!    tem1     3-dimensional array (nz >= 4) to store:
!    tem2     3-dimensional array (nz >= 4) to store:
!    tem3     3-dimensional array (nz >= 4) to store:
!    tem4     3-dimensional array (nz >= 4) to store:
!    tem5     3-dimensional array (nz >= 4) to store:
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

  INTEGER :: nx,ny,nz       ! The number grid points in 3 directions
  INTEGER :: nzsoil            ! The number grid points in the soil

  INCLUDE 'timelvls.inc'

  REAL :: u      (nx,ny,nz) ! Total u-velocity (m/s)
  REAL :: v      (nx,ny,nz) ! Total v-velocity (m/s)
  REAL :: w      (nx,ny,nz) ! Total w-velocity (m/s)

  REAL :: ptprt  (nx,ny,nz) ! Perturbation potential temperature (K)
  REAL :: pprt   (nx,ny,nz) ! Perturbation pressure (Pascal)
  REAL :: qv     (nx,ny,nz) ! Water vapor specific humidity (kg/kg)
!  REAL :: qr     (nx,ny,nz) ! Rain water mixing ratio (kg/kg)

  REAL :: rhostr (nx,ny,nz) ! Base state density rhobar times j3.
  REAL :: ptbar  (nx,ny,nz) ! Base state potential temperature (K)
  REAL :: pbar   (nx,ny,nz) ! Base state pressure (Pascal)
  REAL :: qvbar  (nx,ny,nz) ! Base state water vapor specific
                            ! humidity (kg/kg)

  REAL :: x      (nx)       ! The x-coord. of the physical and
                            ! computational grid.
                            ! Defined at u-point.
  REAL :: y      (ny)       ! The y-coord. of the physical and
                            ! computational grid.
                            ! Defined at v-point.
  REAL :: z      (nz)       ! The z-coord. of the physical and
                            ! computational grid. Defined at w-point.

  REAL :: zp     (nx,ny,nz) ! The height of the terrain.

  REAL :: zpsoil (nx,ny,nzsoil) ! The depth of the soil.

  REAL :: j1     (nx,ny,nz) ! Coordinate transformation
                            ! Jacobian -d(zp)/d(x)
  REAL :: j2     (nx,ny,nz) ! Coordinate transformation
                            ! Jacobian -d(zp)/d(y)
  REAL :: j3     (nx,ny,nz) ! Coordinate transformation

  REAL :: j3soil (nx,ny,nzsoil)   ! Coordinate transformation
                            ! Jacobian  d(zpsoil)/d(zsoil)
  REAL :: j3inv  (nx,ny,nz) ! Inverse of j3
  REAL :: j3soilinv (nx,ny,nzsoil) ! Inverse of j3soil

  REAL :: prcrate(nx,ny)    ! precipitation rate (kg/(m**2*s))

  INTEGER :: nstyps                ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps) ! Soil type at each point
  REAL :: stypfrct(nx,ny,nstyps)! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)        ! Vegetation type at each point
  REAL :: lai    (nx,ny)           ! Leaf Area Index
  REAL :: roufns (nx,ny)           ! Surface roughness
  REAL :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps)  ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps)  ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)  ! Canopy water amount

  REAL :: qvsfc  (nx,ny,0:nstyps)  ! Effective S.H. at sfc.
                                   ! at the ground level (K)
  REAL :: snowdpth(nx,ny)          ! Snow depth

  REAL :: radsw  (nx,ny)    ! Solar radiation flux reaching the surface
  REAL :: rnflx  (nx,ny)    ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)   ! Net solar radiation, SWin - SWout
  REAL :: radlwin(nx,ny)    ! Incoming longwave radiation

  REAL :: cdm    (nx,ny)    ! Drag coefficient for momentum
                            ! defined as scalar point
  REAL :: cdh    (nx,ny)    ! Drag coefficient for heat
  REAL :: cdq    (nx,ny)    ! Drag coefficient for moisture

  REAL :: usflx  (nx,ny)    ! surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx  (nx,ny)    ! surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx (nx,ny)    ! surface flux of heat (K*kg/(m**2*s))
  REAL :: qvsflx (nx,ny)    ! surface flux of moisture (kg/(m**2*s))

  REAL :: pbldpth(nx,ny,nt) ! Planetary boundary layer depth (m)

  REAL :: ptsfc  (nx,ny)    ! Potential temperature at the ground level (K)
  REAL :: psfc   (nx,ny)    ! Surface pressure (Pascal)
  REAL :: wdspsfc(nx,ny)    ! Wind speed at the surface
  REAL :: rhosfc (nx,ny)    ! Base-state air density at the surface
  REAL :: pt1    (nx,ny)    ! Potential temperature at 1st model level AGL
  REAL :: tair   (nx,ny)    ! Temperature of surface air
  REAL :: qvair  (nx,ny)    ! Specific humidity of surface air

  REAL :: tem1   (nx,ny,nz) ! Temporary array
  REAL :: tem2   (nx,ny,nz) ! Temporary array
  REAL :: tem3   (nx,ny,nz) ! Temporary array
  REAL :: tem4   (nx,ny,nz) ! Temporary array
  REAL :: tem5   (nx,ny,nz) ! Temporary array

  REAL :: tsdiffus (nx,ny,nzsoil) ! Temporary array
  REAL :: tem1soil (nx,ny,nzsoil) ! Temporary array
  REAL :: tem2soil (nx,ny,nzsoil) ! Temporary array
  REAL :: tem3soil (nx,ny,nzsoil) ! Temporary array
  REAL :: tem4soil (nx,ny,nzsoil) ! Temporary array


!-----------------------------------------------------------------------
! Local work arrays:
!-----------------------------------------------------------------------

  REAL :: shflx  (nx,ny)  ! Sensible heat flux
  REAL :: lhflx  (nx,ny)  ! Latent heat flux
  REAL :: gflx   (nx,ny)  ! Diffusive heat flux from ground surface to deep soil
  REAL :: ct     (nx,ny)  ! Thermal capacity
  REAL :: evaprg (nx,ny)  ! Evaporation from groud surface
  REAL :: evaprtr(nx,ny)  ! Transpiration of the remaining part (1-delta) of leaves
  REAL :: evaprr (nx,ny)  ! Direct evaporation from the fraction delta
  REAL :: qvsat  (nx,ny)  ! Surface specific humidity at saturation
  REAL :: qvsata (nx,ny)  ! qvsat(tair) (kg/kg)
  REAL :: f34    (nx,ny)  ! f34 and surface resistance

  REAL :: temxy1 (nx,ny)  ! Temporary array
  REAL :: temxy2 (nx,ny)  ! Temporary array
  REAL :: temxy3 (nx,ny)  ! Temporary array
  REAL :: temxy4 (nx,ny)  ! Temporary array

  REAL :: temxpy (nx+ny)  ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, is, nsmth

  REAL :: tema,temb
  REAL :: gumove,gvmove
  REAL :: nsfcstsave
  INTEGER :: tmstrln
  CHARACTER (LEN=256) :: soiloutfl,temchar
  CHARACTER (LEN=80 ) :: timsnd
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'sfcphycst.inc'

  INCLUDE 'mp.inc'            ! Message passing parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF ( sfcphy == 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  If sfcphy = 0, usflx, vsflx, ptsflx and qvsflx are set to zero and
!  will not be used in the mixing subroutines.
!  For safety, we set them to zero.
!
!-----------------------------------------------------------------------
!
    usflx (:,:) = 0.0
    vsflx (:,:) = 0.0
    ptsflx(:,:) = 0.0
    qvsflx(:,:) = 0.0

    RETURN

  END IF

!
!-----------------------------------------------------------------------
!
!  Prepare surface arrays to be used later.
!
!-----------------------------------------------------------------------
!
  IF ( grdtrns == 0 ) THEN
    gumove = 0.0
    gvmove = 0.0
  ELSE
    gumove = umove
    gvmove = vmove
  END IF

  CALL set_acct(sfcphy_acct)

  DO j=1,ny-1
    DO i=1,nx-1

!      roufns(i,j) = 0.03

      tema = 0.25*(u(i+1,j,1)+u(i,j,1)+u(i+1,j,2)+u(i,j,2))+gumove
      temb = 0.25*(v(i,j+1,1)+v(i,j,1)+v(i,j+1,2)+v(i,j,2))+gvmove
      wdspsfc(i,j)=MAX(vsfcmin,                                         &
                   SQRT(tema*tema+temb*temb+w(i,j,2)*w(i,j,2)))

      psfc  (i,j) = 0.5 * ( pprt(i,j,1) + pbar(i,j,1)                   &
                          + pprt(i,j,2) + pbar(i,j,2) )

      rhosfc(i,j) = 0.5 * ( rhostr(i,j,1)*j3inv(i,j,1)                  &
                          + rhostr(i,j,2)*j3inv(i,j,2) )

      pt1(i,j)=ptbar(i,j,2)+ptprt(i,j,2)

      ptsfc (i,j) = tsoil(i,j,1,0)*(p0/psfc(i,j))**rddcp
    END DO
  END DO

  IF ( sfcphy == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Compute the surface fluxes by calling SFCFLX given drag coefficients
!  ground surface potential temperature and equivalent specific humidity.
!
!-----------------------------------------------------------------------
!

    DO j=1,ny-1
      DO i=1,nx-1
        IF ( soiltyp(i,j,1) == 13 .AND.landwtr == 1) THEN
          cdm(i,j) = cdmwtr
          cdh(i,j) = cdhwtr
          cdq(i,j) = cdqwtr
        ELSE
          cdm(i,j) = cdmlnd
          cdh(i,j) = cdhlnd
          cdq(i,j) = cdqlnd
        END IF
      END DO
    END DO

    CALL sfcflx(nx,ny,nz,u,v,w,ptprt,pprt,qv,ptbar,pbar,qvbar,          &
         wdspsfc,rhosfc,ptsfc,qvsfc(1,1,0), cdm,cdh,cdq,                &
         usflx,vsflx,ptsflx,qvsflx)

  ELSE IF ( sfcphy == 2 ) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate stability-dependent drag coefficients by calling DRAGCOEF.
!  Then calculate the fluxes.
!
!-----------------------------------------------------------------------
!

    CALL dragcoef(nx,ny,nz,                                             &
         zp, ptsfc,qvsfc,soiltyp,roufns,wdspsfc,pt1,                    &
         cdm,cdh,cdq)

    CALL sfcflx(nx,ny,nz,u,v,w,ptprt,pprt,qv,ptbar,pbar,qvbar,          &
         wdspsfc,rhosfc,ptsfc,qvsfc(1,1,0), cdm,cdh,cdq,                &
         usflx,vsflx,ptsflx,qvsflx)

  ELSE IF ( sfcphy == 3 .OR. sfcphy == 4 ) THEN
!
!-----------------------------------------------------------------------
!
!  Before calling SFCEBM, we need to prepare some data for the
!  surface model first.
!
!-----------------------------------------------------------------------
!
    DO j=1,ny-1
      DO i=1,nx-1
        tair(i,j) = 0.5 * ( (ptbar(i,j,1)+ptprt(i,j,1))                 &
                          + (ptbar(i,j,2)+ptprt(i,j,2)) )               &
                        * (psfc(i,j)/p0)**rddcp   ! Temperature
        qvair(i,j) = 0.5 * ( qv(i,j,1) + qv(i,j,2) )
      END DO
    END DO

    tsoil(:,:,:,0) = 0.0
    qsoil(:,:,:,0) = 0.0
    wetcanp(:,:,0) = 0.0
    qvsfc(:,:,0)   = 0.0

    usflx (:,:) = 0.0
    vsflx (:,:) = 0.0
    ptsflx(:,:) = 0.0
    qvsflx(:,:) = 0.0

    IF ( soilinitopt == 1 .AND. nstep == 1 ) THEN
      nsfcstsave = nsfcst
      nsfcst = nint(soiltintv/dtsfc)
    END IF

    DO is=nstyp,1,-1
!
!-----------------------------------------------------------------------
!
!  Specified drag coefficients cdm, cdh and cdq.
!
!-----------------------------------------------------------------------
!

      DO j=1,ny-1
        DO i=1,nx-1
          ptsfc (i,j) = tsoil(i,j,1,is)*(p0/psfc(i,j))**rddcp
        END DO
      END DO

      IF ( sfcphy == 3 ) THEN

        DO j=1,ny-1
          DO i=1,nx-1
            IF ( soiltyp(i,j,is) == 13 .AND.landwtr == 1) THEN
              cdm(i,j) = cdmwtr
              cdh(i,j) = cdhwtr
              cdq(i,j) = cdqwtr
            ELSE
              cdm(i,j) = cdmlnd
              cdh(i,j) = cdhlnd
              cdq(i,j) = cdqlnd
            END IF
          END DO
        END DO

      END IF

      IF ( sfcphy == 4 ) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate drag coefficients and surface fluxes for each subtypes
!  within a grid cell. The final fluxes are weighted averaged of all
!  types. Drag coefficients are type dependent and are therefore
!  recalculated for each type.
!
!-----------------------------------------------------------------------
!
        CALL dragcoef(nx,ny,nz,zp,                                      &
             ptsfc,qvsfc(1,1,is),soiltyp(1,1,is),roufns,wdspsfc,pt1,    &
             cdm,cdh,cdq)
      END IF
!
!-----------------------------------------------------------------------
!
!  usflx, vsflx, ptsflx and qvsflx are stored in temxy1, temxy2, temxy3
!  and temxy4, respectively.
!
!-----------------------------------------------------------------------

      CALL sfcflx(nx,ny,nz,u,v,w,ptprt,pprt,qv,ptbar,pbar,qvbar,        &
           wdspsfc,rhosfc,ptsfc,qvsfc(1,1,is), cdm,cdh,cdq,             &
           temxy1,temxy2,temxy3,temxy4)

      DO i=2,nx-1
        DO j=1,ny-1
          usflx(i,j)  = usflx(i,j)  + temxy1(i,j)*stypfrct(i,j,is)
        END DO
      END DO

      DO i=1,nx-1
        DO j=2,ny-1
          vsflx(i,j)  = vsflx(i,j)  + temxy2(i,j)*stypfrct(i,j,is)
        END DO
      END DO

      DO i=1,nx-1
        DO j=1,ny-1
          ptsflx(i,j) = ptsflx(i,j) + temxy3(i,j)*stypfrct(i,j,is)
          qvsflx(i,j) = qvsflx(i,j) + temxy4(i,j)*stypfrct(i,j,is)
        END DO
      END DO

!-----------------------------------------------------------------------
!  Integrate soil model (soilebm) to update tsoil, qvsfc, qsoil,
!  and wetcanp. They will be used by the next step.
!-----------------------------------------------------------------------

      CALL set_acct(soil_acct)

      IF (soilmodel_option == 1) THEN           ! Force-restore scheme

       CALL soilebm(nx,ny,nz,                                            &
            soiltyp(1,1,is),vegtyp,lai,veg,                              &
            tsoil(1,1,1,is),tsoil(1,1,2,is),qsoil(1,1,1,is),             &
            qsoil(1,1,2,is),wetcanp(1,1,is),snowdpth,qvsfc(1,1,is),      &
            wdspsfc,psfc,rhosfc,prcrate,                                 &
            tair,qvair,cdh,cdq,radsw,rnflx,                              &
            shflx,lhflx,gflx,ct,evaprg,evaprtr,evaprr,qvsat,qvsata,f34)

      ELSE IF (soilmodel_option == 2) THEN     ! Implicit scheme

       CALL ousoil(nx,ny,nz,nzsoil,zpsoil,j3soil,j3soilinv,              &
            soiltyp(1,1,is),vegtyp,lai,veg,                              &
            tsoil(1,1,1,is),qsoil(1,1,1,is),                             &
            wetcanp(1,1,is),snowdpth,qvsfc(1,1,is),                      &
            wdspsfc,psfc,rhosfc,prcrate,                                 &
            tair,qvair,cdm,cdh,cdq,radsw,rnflx, radswnet, radlwin,       &
            shflx,lhflx,gflx,ct,evaprg,evaprtr,evaprr,qvsat,qvsata,f34,  &
            tem1soil,tem2soil,tem3soil,tem4soil,tsdiffus,tem1(1,1,1),    &
            tem2(1,1,2),tem2(1,1,3) )

      END IF


!-----------------------------------------------------------------------
!  If we want to use the sensible and latent heat fluxes calculated
!  inside the soil model (soilebm) in the atmospheric model, doing the
!  folloing do loop. Otherwise, the flux calculated in sfcflx will be
!  used, which is based on qvsfc one time step earlier. Otherwise,
!  they should be the same.
!-----------------------------------------------------------------------
!
!     DO i=1,nx-1
!       DO j=1,ny-1
!         ptsflx(i,j) = ptsflx(i,j) + shflx(i,j)*stypfrct(i,j,is)
!         qvsflx(i,j) = qvsflx(i,j) + lhflx(i,j)*stypfrct(i,j,is)
!       END DO
!     END DO

      DO k=1,nzsoil
        DO j=1,ny-1
          DO i=1,nx-1
            tsoil  (i,j,k,0) = tsoil  (i,j,k,0)                         &
                    + tsoil  (i,j,k,is)*stypfrct(i,j,is)
            qsoil  (i,j,k,0) = qsoil  (i,j,k,0)                         &
                    + qsoil  (i,j,k,is)*stypfrct(i,j,is)
          END DO
        END DO
      END DO

      DO j=1,ny-1
        DO i=1,nx-1
          wetcanp(i,j,0) = wetcanp(i,j,0)                               &
                         + wetcanp(i,j,is)*stypfrct(i,j,is)
          qvsfc  (i,j,0) = qvsfc  (i,j,0)                               &
                         + qvsfc  (i,j,is)*stypfrct(i,j,is)
        END DO
      END DO

    END DO  ! loop over is

!-----------------------------------------------------------------------
!  If we want to use the sensible and latent heat fluxes calculated
!  inside the soil model (soilebm) in the atmospheric model, doing the
!  folloing do loop. Otherwise, the flux calculated in sfcflx will be
!  used, which is based on qvsfc one time step earlier. Otherwise,
!  they should be the same.
!
!  Convert sensible and latent heat fluxes stored currently in
!  ptsflx and qvsflx into pt and qv fluxes required by atmospheric model.
!-----------------------------------------------------------------------
!
!   tem = 1.0/lathv
!   DO i=1,nx-1
!     DO j=1,ny-1
!       ptsflx(i,j)=ptsflx(i,j)*(p0/(pbar(i,j,1)+ptbar(i,j,2)          &
!                   +pprt(i,j,1)+pprt(i,j,2))*0.5)**rddcp
!       qvsflx(i,j)=qvsflx(i,j)*tem
!     ENDDO
!   ENDDO

    IF ( soilinitopt == 1 .AND. nstep == 1 ) THEN
      nsfcst = nsfcstsave

      CALL cvttsnd( curtim, timsnd, tmstrln )

      soiloutfl = runname(1:lfnkey)//'.new_soilvar.'//timsnd(1:tmstrln)

      IF( dirname /= ' ' ) THEN
        temchar = soiloutfl
        soiloutfl = dirname(1:ldirnam)//'/'//temchar
      END IF

      IF(mp_opt > 0 .AND. joindmp(FINDX_S) > 0) THEN
        CALL wrtjoinsoil( nx,ny,nzsoil,nstyps, soiloutfl,               &
                 dx,dy,zpsoil,                                          &
                 mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,   &
                 1,1,1,1,1,                                             &
                 tsoil,qsoil,wetcanp,snowdpth,soiltyp)
      ELSE
        CALL wrtsoil( nx,ny,nzsoil,nstyps, soiloutfl,                   &
                 dx,dy,zpsoil,                                          &
                 mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,   &
                 1,1,1,1,1,                                             &
                 tsoil,qsoil,wetcanp,snowdpth,soiltyp)

        CALL soilcntl( nx,ny,nzsoil,zpsoil,soiloutfl,1,1,1,1,1, x,y)

      END IF

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Update the boundary zone fluxes.
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv1dew(usflx,nx,ny,ebc,wbc,1,temxpy)
    CALL mpsendrecv1dns(usflx,nx,ny,nbc,sbc,1,temxpy)

    CALL mpsendrecv1dew(vsflx,nx,ny,ebc,wbc,2,temxpy)
    CALL mpsendrecv1dns(vsflx,nx,ny,nbc,sbc,2,temxpy)

    CALL mpsendrecv1dew(ptsflx,nx,ny,ebc,wbc,0,temxpy)
    CALL mpsendrecv1dns(ptsflx,nx,ny,nbc,sbc,0,temxpy)

    CALL mpsendrecv1dew(qvsflx,nx,ny,ebc,wbc,0,temxpy)
    CALL mpsendrecv1dns(qvsflx,nx,ny,nbc,sbc,0,temxpy)
  END IF

  CALL acct_interrupt(bc_acct)
  CALL bcu2d(nx,ny, usflx,  ebc,wbc,nbc,sbc)
  CALL bcv2d(nx,ny, vsflx,  ebc,wbc,nbc,sbc)
  CALL bcs2d(nx,ny, ptsflx, ebc,wbc,nbc,sbc)
  CALL bcs2d(nx,ny, qvsflx, ebc,wbc,nbc,sbc)
  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
!  Smooth fluxes to remove any 2-dx waves caused by discontinuities
!  in the soil conditions.
!
!-----------------------------------------------------------------------
!
  IF ( smthflx == 1 ) THEN
    DO nsmth=1,numsmth
      CALL smooth9p_nobc( usflx, nx,ny,1,nx,  1,ny-1,temxy1 )
      CALL smooth9p_nobc( vsflx, nx,ny,1,nx-1,1,ny,  temxy1 )
      CALL smooth9p_nobc( ptsflx,nx,ny,1,nx-1,1,ny-1,temxy1 )
      CALL smooth9p_nobc( qvsflx,nx,ny,1,nx-1,1,ny-1,temxy1 )

      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv1dew(usflx,nx,ny,ebc,wbc,1,temxpy)
        CALL mpsendrecv1dns(usflx,nx,ny,nbc,sbc,1,temxpy)

        CALL mpsendrecv1dew(vsflx,nx,ny,ebc,wbc,2,temxpy)
        CALL mpsendrecv1dns(vsflx,nx,ny,nbc,sbc,2,temxpy)

        CALL mpsendrecv1dew(ptsflx,nx,ny,ebc,wbc,0,temxpy)
        CALL mpsendrecv1dns(ptsflx,nx,ny,nbc,sbc,0,temxpy)

        CALL mpsendrecv1dew(qvsflx,nx,ny,ebc,wbc,0,temxpy)
        CALL mpsendrecv1dns(qvsflx,nx,ny,nbc,sbc,0,temxpy)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcu2d(nx,ny, usflx,  ebc,wbc,nbc,sbc)
      CALL bcv2d(nx,ny, vsflx,  ebc,wbc,nbc,sbc)
      CALL bcs2d(nx,ny, ptsflx, ebc,wbc,nbc,sbc)
      CALL bcs2d(nx,ny, qvsflx, ebc,wbc,nbc,sbc)
      CALL acct_stop_inter
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Write all surface related fields for diagnostic purpose.
!
!-----------------------------------------------------------------------
!
  IF ( sfcdiag == 1 ) THEN

    CALL set_acct(output_acct)

    CALL soildiag(nx,ny,nzsoil,x,y,zpsoil,                              &
         soiltyp(1,1,1),vegtyp,lai,roufns,veg,zp(1,1,2),                &
         tsoil(1,1,1,0),qsoil(1,1,1,0),                                 &
         wetcanp(1,1,0),qvsfc(1,1,0),                                   &
         usflx,vsflx,ptsflx,qvsflx,                                     &
         wdspsfc,psfc,rhosfc,prcrate,                                   &
         tair,qvair,cdh,cdq,cdm,                                        &
         radsw,rnflx,                                                   &
         shflx,lhflx,gflx,ct,evaprg,evaprtr,evaprr,qvsat,qvsata,f34,    &
         tem1soil)

  END IF
!
!-----------------------------------------------------------------------
!  Call the subroutine to determine height of PBL top
!-----------------------------------------------------------------------
!
  IF (pbldopt > 0 .AND. tmixopt < 5) THEN

    CALL set_acct(sfcphy_acct)

    CALL pbldepth(nx,ny,nz, u,v,w,ptprt,qv,ptbar, zp,j3, ptsflx,        &
                  pbldpth)

  END IF

  RETURN
END SUBROUTINE sfcphysics
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SFCFLX                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sfcflx(nx,ny,nz,                                             &
           u,v,w,ptprt,pprt,qv,ptbar,pbar,qvbar,                        &
           wdspsfc,rhosfc,ptsfc,qvsfc, cdm,cdh,cdq,                     &
           usflx,vsflx,ptsflx,qvsflx)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the surface momentum, heat and surface moisture fluxes.
!  These fluxes will be passed into turbulent mixing subroutines.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/07/1992.
!
!  MODIFICATION HISTORY:
!
!  2/24/1999 (M. Xue)
!  Simplified codes that calculate the sfc fluxes.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        z component of velocity at a given time level (m/s)
!
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    qvbar    Base state specific humidity (kg/kg)
!
!    wdspsfc  Wind speed at the surface
!    rhosfc   Base-state air density at the surface
!    ptsfc    Potential temperature at the ground level (K)
!    qvsfc    Effective specific humidity at sfc.
!    cdm      Drag coefficient for momentum defined as scalar point
!    cdh      Drag coefficient for heat
!    cdq      Drag coefficient for moisture
!
!  OUTPUT:
!
!    usflx    Surface flux of u-momentum
!    vsflx    Surface flux of v-momentum
!    ptsflx   Surface flux of heat (K*kg/(m**2*s))
!    qvsflx   Surface flux of moisture (K*kg/(m**2*s))
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

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions
  REAL :: u      (nx,ny,nz) ! Total u-velocity (m/s)
  REAL :: v      (nx,ny,nz) ! Total v-velocity (m/s)
  REAL :: w      (nx,ny,nz) ! Total w-velocity (m/s)

  REAL :: ptprt  (nx,ny,nz) ! Perturbation potential temperature (K)
  REAL :: pprt   (nx,ny,nz) ! Perturbation pressure (Pascal)
  REAL :: qv     (nx,ny,nz) ! Water vapor specific humidity (kg/kg)

  REAL :: ptbar  (nx,ny,nz) ! Base state potential temperature (K)
  REAL :: pbar   (nx,ny,nz) ! Base state pressure (Pascal)
  REAL :: qvbar  (nx,ny,nz) ! Base state specific humidity (kg/kg)

  REAL :: ptsfc  (nx,ny)    ! Potential temperature at the ground level (K)
  REAL :: rhosfc (nx,ny)    ! Surface air density (kg/m**3)
  REAL :: wdspsfc(nx,ny)    ! Wind speed at the surface
  REAL :: qvsfc  (nx,ny)    ! Effective specific humidity at surface

  REAL :: cdm    (nx,ny)    ! Drag coefficient for surface momentum flux
  REAL :: cdh    (nx,ny)    ! Drag coefficient for surface heat flux
  REAL :: cdq    (nx,ny)    ! Drag coefficient for surface moisture flux

  REAL :: usflx  (nx,ny)    ! surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx  (nx,ny)    ! surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx (nx,ny)    ! surface flux of heat (K*kg/(m**2*s))
  REAL :: qvsflx (nx,ny)    ! surface flux of moisture (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: usuf, vsuf,qvsfc0
  REAL :: tema
  REAL :: gumove,gvmove
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!  Include file globcst.inc passes cdmlnd,cdmwtr, cdhlnd,cdmwtr,
!  cdqlnd,cdmwtr and UMOVE and VMOVE into this subroutine.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
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
!
!-----------------------------------------------------------------------
!
!  Compute the surface flux of u-momentum, i.e. the surface
!  drag in x-direction, according to a simple drag law:
!
!  usflx = rhobar *bar( (w'u') ) = rhobar *Cdm*V*u
!
!  where Cdm is the drag coefficient, V the wind speed at the first
!  scalar level above ground, and u the x-component of wind at the
!  same level.
!
!  usflx is defined one-half dz below the u-point, i.e. at ground level.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=2,nx-1
      usuf = 0.5*(u(i,j,2)+u(i,j,1))+gumove
      usflx(i,j)=0.125*(cdm(i,j)+cdm(i-1,j))*                           &
                 (rhosfc(i,j)+rhosfc(i-1,j))*                           &
                 (wdspsfc(i,j)+wdspsfc(i-1,j))*                         &
                 sign(1.0, usuf)*MAX(abs(usuf),vsfcmin)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Compute the surface flux of v-momentum, i.e. the surface
!  drag in y-direction, according to a simple drag law:
!
!  vsflx = - rhobar *bar( (w'v') ) = rhobar *Cdm*V*v
!
!  where Cdm is the drag coefficient, V the wind speed at the first
!  scalar level above ground, and v the y-component of wind at the
!  same level.
!
!  vsflx is defined one-half dz below the v-point, i.e. at ground level.
!
!-----------------------------------------------------------------------
!
  DO j=2,ny-1
    DO i=1,nx-1
      vsuf = gvmove + 0.5*(v(i,j,2)+v(i,j,1))
      vsflx(i,j)=0.125*(cdm(i,j)+cdm(i,j-1))*                           &
                 (rhosfc(i,j)+rhosfc(i,j-1))*                           &
                 (wdspsfc(i,j)+wdspsfc(i,j-1))*                         &
                 sign(1.0,vsuf)*MAX(abs(vsuf),vsfcmin)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Compute the surface heat flux according to a simple drag law:
!
!  ptsflx = - rhobar *bar( (w' pt') ) = rhobar *Cdh*V*(pt-ptsfc0)
!
!  where Cdh is the drag coefficient, V the wind speed and pt the
!  potential temperatureat at the first scalar level above ground.
!
!  ptsfc0 is defined as the potential temperature at ground leval.
!  ptsfc0 = ptbar at surface
!
!  ptsflx is defined one-half dz below the scalar point, i.e. at
!  ground level.
!
!  Similarly, the vertical moisture flux at the surface is defined as
!
!  qvsflx = - rhobar *bar( (w' qv') ) = rhobar *Cdq*V*(qv-qvsfc0)
!
!  where qv is the specific humidity at the first model scalar level.
!  qvsfc0 is defined as the mixing ratio at the ground and taken as
!  the base state moisture, qvsfc0 = qvbar, at ground level.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
      IF ( sfcphy == 1 .OR. sfcphy == 3 ) THEN
        qvsfc0=0.5*(qvbar(i,j,1)+qvbar(i,j,2))
      ELSE
        qvsfc0=qvsfc(i,j)
      END IF
      tema = rhosfc(i,j)*MAX(wdspsfc(i,j),vsfcmin)
      ptsflx(i,j)=cdh(i,j)*tema*(ptprt(i,j,2)+ptbar(i,j,2)-ptsfc(i,j))
      qvsflx(i,j)=cdq(i,j)*tema*(qv(i,j,2)-qvsfc0)
    END DO
  END DO

  RETURN
END SUBROUTINE sfcflx
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE DRAGCOEF                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dragcoef(nx,ny,nz,                                           &
           zp, ptsfc,qvsfc,soiltyp,roufns,wdspsfc,pt1,                  &
           cdm,cdh,cdq)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the surface momentum, heat and surface moisture fluxes
!  using a stability dependent formulation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/07/1992.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  6/4/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  9/14/1992 (M. Xue)
!  Different drag coefficients defined for momentum, temperature and
!  moisture. Definition of ptsfc and qvsfc changed. qvbar added to
!  the argument list.
!
!  9/28/1992 (M. Xue)
!  Sign correction in all the flux formulations.
!
!  8/28/1993 (M. Xue and Hao jin)
!  Modified the flux formulations considering the terrain effect.
!
!  9/10/1993 (V. Wong, X. Song and N. Lin)
!  Stability dependent formulation used for usflx, vsflx, and ptsflx.
!  The different flux formulations come from subroutines USTARC and
!  PTSTARC.
!
!  9/04/1994 (K. Brewster, V. Wong and X. Song)
!  Facelift.
!
!  12/7/1994 (Y. Liu)
!  Cleaned up an empty DO loop with a label 360.
!
!  01/28/1995 (V. Wong and X. Song)
!  Add the option of stability and roughness dependent surface layer
!  parameterization for both land and sea surfaces.
!
!  03/02/1995 (Y. Liu)
!  Added the minimum surface total wind speed to guarantee the basic
!  heat and moisture fluxes.
!
!  02/07/96 (V.Wong and X.Song)
!  Changed the calculation of roughness zo over the sea according to
!  Clarke(1970).  The new formulation is insensitive to the value of
!  the depth. of the first model layer above the ground.
!  Added kvwtr to denote the Von Karman constant over the sea.
!  Set a lower limiter, zolimit, for zo, and an upper limiter, z1limit,
!  for depth of the surface layer z1.
!
!  03/17/1997 (Yuhe Liu, Vince Wong and Jing Tian)
!  Added an option for constant cdh (cdq) over water.
!
!  05/29/97 (V. Wong and X. Tan)
!  Modified the formulation considering the height of the surface
!  layer z1 may equal zero.
!
!  02/25/1999 (M. Xue)
!  Isolated from original SFCFLXSD routine. The above listing is for
!  SFCFLXSD.
!
!  01/09/2002 (Y. Wang)
!  Changed the DO loop lower bound from i=2, j=2 to i=1, j=1 in
!  the calculation of drag coefficient for momentum.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!
!    ptsfc    Potential temperature at the ground level (K)
!    qvsfc    Effective specific humidity at sfc
!    soiltyp  Soil type at each point
!    roufns   Surface roughness
!
!    wdspsfc  Surface wind speed (at first level about ground) (m/s)
!    pt1      Potential temperature at the first model level AGL
!
!  OUTPUT:
!
!    cdm      Drag coefficient for surface momentum flux (nondimensional)
!    cdh      Drag coefficient for surface heat flux (nondimensional)
!    cdq      Drag coefficient for surface moisture flux (nondimensional)
!
!  Local work arrays:
!
!    c_u      Cu
!    c_pt     Cpt
!    c_uneu   Cu for neutral stratification
!    c_ptneu  Cpt for neutral stratification
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of staggered grid.
  REAL :: ptsfc (nx,ny)        ! Potential temperature at the ground level (K)
  REAL :: qvsfc (nx,ny)        ! Effective surface specific humidity

  INTEGER :: soiltyp(nx,ny)    ! Soil type at each point
  REAL :: roufns(nx,ny)        ! Surface roughness length

  REAL :: wdspsfc(nx,ny)       ! Wind speed at the surface
  REAL :: pt1   (nx,ny)        ! Potential temperature at the first model level AGL

  REAL :: cdm   (nx,ny)        ! Drag coefficient for surface momentum flux
  REAL :: cdh   (nx,ny)        ! Drag coefficient for surface heat flux
  REAL :: cdq   (nx,ny)        ! Drag coefficient for surface moisture flux

  REAL :: c_u   (nx,ny)        ! Cu
  REAL :: c_pt  (nx,ny)        ! Cpt
  REAL :: c_uneu(nx,ny)        ! Cu for neutral stratification
  REAL :: c_ptneu(nx,ny)       ! Cpt for neutral stratification
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: check
  REAL :: z1,z1drou,z1droup
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!  Include file globcst.inc passes cdq (moisture transfer coeffcient)
!  and UMOVE and VMOVE into this subroutine.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny-1
    DO i=1,nx-1
!
! Calculate C_u and C_pt at the neutral state
!
      z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
      z1  = MIN(z1,z1limit)
      z1drou = z1/roufns(i,j)

      z1droup = (z1+roufns(i,j))/roufns(i,j)

      c_uneu (i,j) = kv/LOG(z1droup)
      c_ptneu(i,j) = c_uneu(i,j) /prantl0l


    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  A general procedure for computing the drag coefficients:
!
!  1.  Obtain the following quantities:
!      a. The first scalar level above ground level (z1)
!      b. Potential temperature at the height of surface roughness
!         (roufns) - given.
!      c. Potential temperature (pt1) at the first scalar level above
!         ground.
!      d. C_u and C_pt at neutral state in the land case from
!             c_uneu =kv/log10(z1/roughness)
!             c_ptneu=c_uneu/prandtl_number
!         In the water case, call subroutines CUNEUWTR and CPTNEUWTR.
!
!  2.  Given the information from (1), compute the Bulk Richardson
!      number, defined as
!
!         BulkRi = (g/ptsfc) * ( (pt1-ptsfc)*(z1-roufns) / V**2 )
!      where V is the wind speed at the surface.
!
!  3.  After calculation of the BulkRi, then:
!      a. Given the BulkRi and the information in (1), compute C_u at
!         both land and water cases (see subroutines CUC and CUCWTR).
!
!      b. From C_u and C_pt, compute draf coefficients cdh and cdq.
!
!-----------------------------------------------------------------------
!
!  Calculate drag coefficient (defined at scalar point) for momentum
!  fluxes.
!
!-----------------------------------------------------------------------
!
  CALL cuc(nx,ny,nz,1,nx-1,1,ny-1,zp,                                   &
           roufns,wdspsfc,ptsfc,pt1,c_uneu, c_u)

  IF (wtrexist == 1) THEN

    CALL cuneuwtr(nx,ny,nz,1,nx-1,1,ny-1,soiltyp,wdspsfc,               &
                  c_uneu)

    CALL cucwtr  (nx,ny,nz,1,nx-1,1,ny-1,zp,soiltyp,wdspsfc,            &
                  ptsfc,pt1,c_uneu, c_u)

  END IF

  DO j=1,ny-1
    DO i=1,nx-1
      check=ptsfc(i,j)-pt1(i,j)
!
!    Put lower and upper limits on C_u
!
      IF(check >= 0.0) THEN  ! Unstable case
        c_u(i,j)=MIN( c_u(i,j), 2.0*c_uneu(i,j))
      ELSE                   ! Stable case
        c_u(i,j)=MAX( c_u(i,j), blimit*c_uneu(i,j))
      END IF
!
!-----------------------------------------------------------------------
!
!  c_dm=c_U**2
!
!  Change Von Karman constant from 0.4 to 0.35 and since the drag
!  coefficient Cdm is proportional to the square of the constant,
!  the final usflx should be multiplied by (0.35/0.40)**2 = 0.765625
!
!-----------------------------------------------------------------------
!
      cdm(i,j) = 0.765625 * c_u(i,j) * c_u(i,j)

    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  From C_u and C_pt, compute C_dh = C_u * C_pt and pt flux,
!
!-----------------------------------------------------------------------
!
  CALL cptc(nx,ny,nz,1,nx-1,1,ny-1,zp,roufns,wdspsfc,                   &
            ptsfc,pt1,c_ptneu,c_pt,c_u)

  IF ( wtrexist == 1 .AND. cdhwtropt == 0 ) THEN

    CALL cptneuwtr(nx,ny,nz,1,nx-1,1,ny-1,soiltyp,wdspsfc,c_uneu,       &
                   c_ptneu)

    CALL cptcwtr(nx,ny,nz,1,nx-1,1,ny-1,zp,soiltyp,wdspsfc,             &
                 ptsfc,pt1,c_ptneu, c_pt)

  END IF

  DO j=1,ny-1
    DO i=1,nx-1
      check=ptsfc(i,j)-pt1(i,j)
!
!    Impose lower and upper limits on C_u and C_pt.
!

!      IF (soilmodel_option == 1) THEN

      IF ( soiltyp(i,j) == 13 .AND. cdhwtropt == 1 ) THEN
        cdh(i,j) = cdhwtr
      ELSE
        IF(check >= 0) THEN
          c_pt(i,j) = MIN(c_pt(i,j),3.3333*c_ptneu(i,j) )
        ELSE
          c_pt(i,j) = MAX(c_pt(i,j),blimit*c_ptneu(i,j) )
        END IF

        cdh(i,j) = c_u(i,j)*c_pt(i,j)
      END IF
!
!-----------------------------------------------------------------------
!
!  Change Von Karman constant from 0.4 to 0.35 and since the drag
!  coefficient Cdh is proportional to the square of the constant,
!  the final Cdh should be multiplied by (0.35/0.40)**2 = 0.765625
!
!-----------------------------------------------------------------------
!
      cdh(i,j) = 0.765625*cdh(i,j)
      cdq(i,j) = 0.7*cdh(i,j)

!      ELSE IF (soilmodel_option == 2) THEN
!        IF ( soiltyp(i,j) == 13 .AND. cdhwtropt == 1 ) THEN
!          cdh(i,j) = cdhwtr
!        ELSE
!          cdh(i,j) = c_pt(i,j)
!          cdq(i,j) = 0.7*cdh(i,j)
!        END IF
!      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE dragcoef
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    SUBROUTINE CUC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cuc(nx,ny,nz,ibgn,iend,jbgn,jend,zp,roufns,wspd,ptsfc,       &
           pt1, c_uneu,c_u)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute C_u (friction velocity / wind speed). The quantity C_U is used by
!  the subroutine SFCFLX to obtain surface fluxes for both unstable
!  and stable cases.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong, X. Song and N. Lin
!  9/10/1993
!
!  For the unstable case (Bulk Richardson number bulkri < 0 ), the
!  formulation is based on the paper "On the Analytical Solutions of
!  Flux-Profile relationships for the Atmospheric Surface Layer" by
!  D.W. Byun in J. of Applied Meteorlogy, July 1990, pp. 652-657.
!
!  For the stable case, the formulation is based on the
!  paper "A Short History of the Operational PBL - Parameterization
!  at ECMWF" by J.F.Louis, M. Tiedtke and J.F. Geleyn in " Workshop
!  on Planetary Boundary Layer Parameterization", a publication
!  by the European Centre for Medium Range Weather Forecasts,
!  25-27 Nov. 1981.
!
!  MODIFICATION HISTORY:
!
!  9/04/1994 (K. Brewster, V. Wong and X. Song)
!  Facelift.
!
!  2/27/95 (V. Wong and X. Song)
!
!  02/07/96 (V.Wong and X.Song)
!  Set an upper limiter, z1limit, for depth of the surface layer z1.
!
!  05/01/97 (V. Wong and X. Tan)
!  Changed the computation of stabp, use the formulation similar to
!  one used by Bynn(1990), the momentum and thermal roughness
!  lengths are different.
!
!  05/29/97 (V. Wong and X. Tan)
!  Modified the formulation considering the height of the surface
!  layer z1 may equal zero.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins.
!    iend     i-index where evaluation ends.
!    jbgn     j-index where evaluation begins.
!    jend     j-index where evaluation ends.
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!    roufns   Surface roughness
!
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!    ptsfc    Potential temperature at the ground level (K)
!    pt1      Potential temperature at the 1st scalar point above
!             ground level, (K)
!
!    c_uneu   Friction velocity at neutral state
!
!  OUTPUT:
!
!    ustar    Friction velocity
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

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins.
  INTEGER :: iend              ! i-index where evaluation ends.
  INTEGER :: jbgn              ! j-index where evaluation begins.
  INTEGER :: jend              ! j-index where evaluation ends.

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate
                               ! defined at w-point of staggered grid.
  REAL :: roufns(nx,ny)        ! Surface roughness length

  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL :: ptsfc(nx,ny)         ! Potential temperature at the ground
                               ! level (K)
  REAL :: pt1   (nx,ny)        ! Potential temperature at the
                               ! 1st scalar
                               ! point above ground level, (K)
  REAL :: c_uneu(nx,ny)      ! Frictional velocity (m/s) at neutral

  REAL :: c_u(nx,ny)        ! Frictional velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j

  REAL :: z1                   ! The height of 1st scalar point above
                               ! ground level (m)
  REAL :: dz0                  ! z1-roufns, where roufns is defined as
                               ! surface roughness length
  REAL :: bulkri               ! Bulk Richardson number

  REAL :: stabp                ! Monin-Obukhov STABility Parameter
                               ! (zeta)
  REAL :: x1,x0,psim           ! Intermediate parameters needed
  REAL :: z1drou,qb3pb2
  REAL :: c7,c8
  REAL :: sb,qb,pb,thetab,tb   ! During  computations
  REAL :: a,b,c,d
  REAL :: tempan,sqrtqb

  REAL :: zt                   ! Thermal roughness length
  REAL :: dzt
  REAL :: ztdrou

  REAL :: z1droup, ztdroup

  REAL, PARAMETER :: epsilon = 1.0E-6

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
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
!  Following Byun (1990).
!
!-----------------------------------------------------------------------
!
  DO j=jbgn,jend
    DO i=ibgn,iend

      z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
      z1  = MIN(z1,z1limit)

      zt = ztz0*roufns(i,j)           ! Original Zt formulation

!---------------------------------------------------------------------
!     Test of new thermal roughness equation (Chen/Dudhia 2001) (JAB)
!---------------------------------------------------------------------
! The following was commented out temporarily because variable psim
! is used before it is initilized. It need to be taken care of later.
!
!      dz0 = z1-roufns(i,j)
!      z1droup = (z1+roufns(i,j))/roufns(i,j)
!
!      bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dz0/                &
!               (wspd(i,j)*wspd(i,j))
!
!      IF (bulkri <= 0.0) THEN
!        c_u(i,j) =kv/(LOG(z1droup)-psim)
!      ELSE
!        a=kv/LOG(z1droup)
!        b=5.0
!        d=5.0
!        c=SQRT(1.0+d*bulkri)
!        c_u(i,j) = a/SQRT(1.0+2.0*b*bulkri/c)
!      END IF
!
!       zt = (0.4*c_u(i,j)/0.000024) + 100.0
!
!       IF (abs(zt) > epsilon ) THEN
!         zt = 1.0/zt
!       ELSE IF (abs(zt) <= epsilon ) THEN
!         zt = ztz0*roufns(i,j)
!       END IF
!
!----------------------------------------------------------------------

      dzt = z1-zt
      ztdrou = z1/zt
      ztdroup = (z1+zt)/zt

      dz0 = z1-roufns(i,j)
      z1drou = z1/roufns(i,j)
      z1droup = (z1+roufns(i,j))/roufns(i,j)

      bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dz0/                &
               (wspd(i,j)*wspd(i,j))

      IF (bulkri <= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  Unstable case: See equations (28)-(34) in Byun (1990).
!
!-----------------------------------------------------------------------
!
        bulkri = MAX (bulkri,-10.0)

        sb =bulkri/prantl0l

        qb=oned9*(c1l+c2l*sb*sb)
        pb=oned54*(c3l+c4l*sb*sb)

        qb3pb2=qb**3-pb*pb
        c7 = (z1*dzt*LOG(z1droup)*LOG(z1droup))/(dz0*dz0*LOG(ztdroup))

        IF( qb3pb2 >= 0.0 ) THEN

          sqrtqb = SQRT(qb)
          tempan = MAX( -1.0, MIN( 1.0, pb/(sqrtqb**3) ) )

          thetab=ACOS(tempan)
          stabp =c7*(-2.0*sqrtqb*COS(thetab/3.0)+c5l)

        ELSE

          tb    =(SQRT(-qb3pb2)+ABS(pb))**oned3
          stabp =c7*(-(tb+qb/tb)+c5l)

        END IF
!
!-----------------------------------------------------------------------
!
!  According to equation (14) in Byun (1990).
!
!-----------------------------------------------------------------------
!
        c8=gammaml*stabp
        x1=(1. - c8)**0.25
        x0=(1. - c8/z1drou)**0.25

        psim=2.0*LOG((1.0+x1)/(1.0+x0))+LOG((1.+x1*x1)/(1.+x0*x0))-     &
             2.0*ATAN(x1)+2.0*ATAN(x0)

!
!-----------------------------------------------------------------------
!
!  Compute C_u via equation (10) in Byun (1981).
!
!-----------------------------------------------------------------------
!
        c_u(i,j) =kv/(LOG(z1droup)-psim)

      ELSE
!
!-----------------------------------------------------------------------
!
!  Stable case:
!
!-----------------------------------------------------------------------
!
        a=kv/LOG(z1droup)
        b=5.0
        d=5.0
        c=SQRT(1.0+d*bulkri)

        c_u(i,j) = a/SQRT(1.0+2.0*b*bulkri/c)

      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cuc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   SUBROUTINE CPTC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cptc(nx,ny,nz,ibgn,iend,jbgn,jend,                           &
           zp,roufns,wspd,ptsfc,pt1, c_ptneu,c_pt,c_u)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute c_pt (a nondimensional temperature scale)
!
!  The quantity c_pt is used by the subroutine SFCFLX to obtain
!  surface fluxes for both the unstable and stable cases.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong, X. Song and N. Lin
!  9/10/1993
!
!  For the unstable case (Bulk Richardson number bulkri < 0 ), the
!  formulation is based on the paper "On the Analytical Solutions of
!  Flux-Profile relationships for the Atmospheric Surface Layer" by
!  D.W. Byun in J. of Applied Meteorlogy, July 1990, pp. 652-657.
!
!  For stable case, the formulation is based on the paper "A Short
!  History of the Operational PBL - Parameterization at ECMWF" by
!  J.F.Louis, M. Tiedtke and J.F. Geleyn in "Workshop on Planetary
!  Boundary Layer Parameterization", a publication by European Centre
!  for Medium Range Weather Forecasts, 25-27 Nov. 1981.
!
!  MODIFICATION HISTORY:
!
!  9/04/1994 (K. Brewster, V. Wong and X. Song)
!  Facelift.
!
!  2/27/95 (V. Wong and X. Song)
!
!  02/07/96 (V.Wong and X.Song)
!  Set an upper limiter, z1limit, for depth of the surface layer z1.
!
!  05/01/97 (V. Wong and X. Tan)
!  Changed the computation of temperature scale
!
!  05/29/97 (V. Wong and X. Tan)
!  Modified the formulation considering the height of the surface
!  layer z1 may equal zero.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins.
!    iend     i-index where evaluation ends.
!    jbgn     j-index where evaluation begins.
!    jend     j-index where evaluation ends.
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!    roufns   Surface roughness
!
!    ptsfc    Potential temperature at the ground level (K)
!    pt1      Potential temperature at the 1st scalar point above
!             ground level, (K)
!
!    c_ptneu  Temperature scale (K) at neutral state, defined by
!             surface heat flux / friction velocity
!
!  OUTPUT:
!
!    c_pt     Temperature scale (K), defined by
!             surface heat flux / friction velocity
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

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins
  INTEGER :: iend              ! i-index where evaluation ends
  INTEGER :: jbgn              ! j-index where evaluation begins
  INTEGER :: jend              ! j-index where evaluation ends

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate
                               ! defined at w-point of staggered grid.
  REAL :: roufns(nx,ny)        ! Surface roughness length

  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL :: ptsfc(nx,ny)         ! Potential temperature at the
                               ! ground level (K)
  REAL :: pt1   (nx,ny)        ! Potential temperature at the
                               ! 1st scalar
                               ! point above ground level, (K)
  REAL :: c_ptneu(nx,ny)       ! Normalized temperature scale (K)
                               ! at neutral state

  REAL :: c_pt  (nx,ny)        ! Temperature scale (K)

  REAL :: c_u   (nx,ny)        ! Friction velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j

  REAL :: z1                  ! The height of 1st scalar point above
                              ! ground level (m)
  REAL :: dz0                 ! z1-roufns, where roufns is defined as
                              ! surface roughness length
  REAL :: bulkri              ! Bulk Richardson number

  REAL :: stabp               ! Monin-Obukhov STABility Parameter
                              ! (zeta)

  REAL :: y1,y0,psih          ! Intermediate parameters needed
  REAL :: z1drou,qb3pb2
  REAL :: c7,c8
  REAL :: sb,qb,pb,thetab,tb  ! During  computations
  REAL :: a,b,c,d
  REAL :: tempan,sqrtqb

  REAL :: zt                  ! Thermal roughness length
  REAL :: dzt,ztdrou

  REAL :: z1droup,ztdroup

  REAL, PARAMETER :: epsilon = 1.0E-6

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
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
!  Following Byun (1990).
!
!-----------------------------------------------------------------------
!
  DO j=jbgn,jend
    DO i=ibgn,iend

      z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
      z1  = MIN(z1,z1limit)

!      zt = ztz0*roufns(i,j)

!---------------------------------------------------------------------
!     Test of new thermal roughness equation (Chen/Dudhia 2001) (JAB)
!---------------------------------------------------------------------
       zt = (0.4*c_u(i,j)/0.000024) + 100.0

       IF (abs(zt) > epsilon ) THEN
         zt = 1.0/zt
       ELSE IF (abs(zt) <= epsilon ) THEN
         zt = ztz0*roufns(i,j)
       END IF

!----------------------------------------------------------------------

      dzt = z1-zt
      ztdrou = z1/zt
      ztdroup = (z1+zt)/zt

      dz0 = z1-roufns(i,j)
      z1drou = z1/roufns(i,j)
      z1droup = (z1+roufns(i,j))/roufns(i,j)

      bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dz0/                &
               (wspd(i,j)*wspd(i,j))


!-------------------------------------------------------------------
!     Soilmodel_option = 2 (Implicit)
!--------------------------------------------------------------------

!     IF (soilmodel_option == 2) THEN
!
!     IF (bulkri <= 0.0) THEN
!---------------------------------------------------------------------
!     Unstable case
!---------------------------------------------------------------------
!
!        stabp = 1.0 - ( (15.0 * bulkri) /                         &
!          ( 1.0 + 7.5 * (10.0*kv*kv/(LOG(z1drou)*LOG(ztdrou))) *  &
!          sqrt(-bulkri * z1drou) ) )
!
!     ELSE
!---------------------------------------------------------------------
!     Stable case
!---------------------------------------------------------------------
!
!        stabp = EXP(-bulkri)
!
!     END IF
!
!     c_pt(i,j) = kv*kv*stabp / (LOG(z1drou)*LOG(ztdrou))
!
!     END IF      !soilmodel_option = 2
!-------------------------------------------------------------------



!-------------------------------------------------------------------
!      Soilmodel_option = 1  (Force-Restore)
!------------------------------------------------------------------

!      IF (soilmodel_option == 1) THEN


      IF (bulkri <= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  Unstable case: See equations (28)-(34) in Byun (1990).
!
!-----------------------------------------------------------------------
!
        bulkri = MAX (bulkri,-10.0)

        sb =bulkri/prantl0l

        qb=oned9*(c1l+c2l*sb*sb)
        pb=oned54*(c3l+c4l*sb*sb)

        qb3pb2=qb**3-pb*pb
        c7 = (z1*dzt*LOG(z1droup)*LOG(z1droup))/(dz0*dz0*LOG(ztdroup))

        IF( qb3pb2 >= 0.0 ) THEN

          sqrtqb = SQRT(qb)
          tempan = MAX( -1.0, MIN( 1.0, pb/(sqrtqb**3) ) )

          thetab=ACOS(tempan)
          stabp =c7*(-2.0*SQRT(qb)*COS(thetab/3.0)+c5l)

        ELSE

          tb    =(SQRT(-qb3pb2)+ABS(pb))**oned3
          stabp =c7*(-(tb+qb/tb)+c5l)

        END IF
!
!-----------------------------------------------------------------------
!
!  According to equation (15) in Byun (1990).
!
!-----------------------------------------------------------------------
!
        c8=gammahl*stabp
        y1=SQRT(1.0 - c8)
        y0=SQRT(1.0 - c8/ztdrou)

        psih=2.0*LOG((y1+1.0)/(y0+1.0))
!
!-----------------------------------------------------------------------
!
!  Compute c_pt via equation (11) in Byun (1981).
!
!-----------------------------------------------------------------------
!
        c_pt(i,j)=kv / (prantl0l*(LOG(ztdroup)-psih))

      ELSE
!
!-----------------------------------------------------------------------
!
!  Stable case: See Louis et al (1981).
!
!-----------------------------------------------------------------------
!
        a=kv/LOG(ztdroup)
        b=5.0
        d=5.0
        c=SQRT(1.0+d*bulkri)

        c_pt(i,j) = SQRT(1.0+2.0*b*bulkri/c)
        c_pt(i,j) = a*c_pt(i,j)/(prantl0l*(1.0+3.0*b*bulkri*c))

      END IF
!    END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cptc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CUCWTR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cucwtr(nx,ny,nz,ibgn,iend,jbgn,jend,zp,soiltyp,              &
           wspd,ptsfc,pt1,c_uwtrneu,c_uwtr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute C_u (friction velocity / wind speed) at sea case. The quantity
!  C_uwtr is used by the subroutine SFCFLX to obtain surface fluxes for both
!  unstable and stable cases.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong and X. Song
!  8/04/1994
!
!  MODIFICATION HISTORY:
!
!  2/27/95 (V.W. and X.S.)
!
!  1/12/96 (V.W. and X.S.)
!  Changed the calculation related to zo.
!  Added kvwtr to denote the Von Karman constant over the sea;
!  Set a lower limiter for zo, zolimit, and an upper limiter for z1, z1limit.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins.
!    iend     i-index where evaluation ends.
!    jbgn     j-index where evaluation begins.
!    jend     j-index where evaluation ends.
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!    ptsfc    Potential temperature at the ground level (K)
!    pt1      Potential temperature at the 1st scalar point above
!             ground level, (K)
!
!  OUTPUT:
!
!    c_uwtr   Friction velocity
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

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins.
  INTEGER :: iend              ! i-index where evaluation ends.
  INTEGER :: jbgn              ! j-index where evaluation begins.
  INTEGER :: jend              ! j-index where evaluation ends.

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of staggered grid.
  INTEGER :: soiltyp(nx,ny)    ! Soil type at each point.

  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL :: ptsfc(nx,ny)         ! Potential temperature at the ground level
                               ! (K)
  REAL :: pt1   (nx,ny)        ! Potential temperature at the 1st scalar
                               ! point above ground level, (K)
  REAL :: c_uwtrneu(nx,ny)        ! Frictional velocity (m/s)
  REAL :: c_uwtr(nx,ny)        ! Frictional velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: z1                   ! The height of 1st scalar point above
                               ! ground level (m)
  REAL :: zo                   ! Defined as surface momentum roughness
                               ! length
  REAL :: zt                   ! Defined as surface heat transfer
                               ! roughness length
  REAL :: dzo                  ! z1-zo
  REAL :: dzt                  ! z1-zt
  REAL :: z1dzo                ! z1/zo
  REAL :: z1dzt                ! z1/zt
  REAL :: xcdn                 ! Cdn (sea)
  REAL :: xcdh                 ! Hot-wired value of Cdh (sea)
  REAL :: bulkri               ! Bulk Richardson number

  REAL :: stabp                ! Monin-Obukhov STABility Parameter
                               ! (zeta)

  REAL :: x1,psim           ! Intermediate parameters needed
  REAL :: qb3pb2
  REAL :: c7,c8
  REAL :: sb,qb,pb,thetab,tb   ! during  computations
  REAL :: a,b,c,d
  REAL :: tempan,sqrtqb
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  xcdh = 1.1E-3

  DO j=jbgn,jend
    DO i=ibgn,iend

      IF (soiltyp (i,j) == 13) THEN

        xcdn = (0.4+0.079*wspd(i,j)) * 1.e-3
        z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
        z1  = MIN(z1,z1limit)
!
! calculate zo and zt
!
        zo =0.032*xcdn*wspd(i,j)*wspd(i,j)/9.8
        zo = MAX(zo,zolimit)
        zt = z1 * EXP( -kvwtr*SQRT(xcdn)/(prantl0w*xcdh) )

        dzo = z1-zo
        z1dzo = z1/zo
        dzt = z1-zt
        z1dzt = z1/zt

        bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dzo*dzo/          &
                 (dzt*wspd(i,j)*wspd(i,j))

        IF (bulkri <= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  Unstable case: A modified formulation, which is similar to
!  equations (28)-(34) in Byun (1990).
!
!-----------------------------------------------------------------------
!
          bulkri = MAX (bulkri,-10.0)

          sb =bulkri/prantl0w

          qb=oned9*(c1w+c2w*sb*sb)
          pb=oned54*(c3w+c4w*sb*sb)

          qb3pb2=qb**3-pb*pb
          c7 = (z1*dzt/(dzo*dzo))*(ALOG(z1dzo)*ALOG(z1dzo)/ALOG(z1dzt))

          IF( qb3pb2 >= 0.0 ) THEN

            sqrtqb = SQRT(qb)
            tempan = MAX( -1.0, MIN( 1.0, pb/(sqrtqb**3) ) )

            thetab=ACOS(tempan)
            stabp =c7*(-2.0*SQRT(qb)*COS(thetab/3.0)+c5w)

          ELSE

            tb    =(SQRT(-qb3pb2)+ABS(pb))**oned3
            stabp =c7*(-(tb+qb/tb)+c5w)

          END IF
!
!-----------------------------------------------------------------------
!
!  According to a modified equation, which is similar to equation (14)
!  in Byun (1990).
!
!-----------------------------------------------------------------------
!
          c8=gammamw*stabp
          x1=(1. - c8)**0.25

          psim=2.0*ALOG(0.5*(1.0+x1))+ALOG(0.5*(1.+x1*x1))-             &
               2.0*ATAN(x1)+ASIN(1.)

!
!-----------------------------------------------------------------------
!
!  Compute c_uwtr via equation (10) in Byun (1981).
!
!-----------------------------------------------------------------------
!
          c_uwtr(i,j) =kv/(LOG(z1dzo)-psim)

        ELSE
!
!-----------------------------------------------------------------------
!
!  Stable case: See Louis et al (1981).
!
!-----------------------------------------------------------------------
!
          a=kv/LOG(z1dzo)
          b=5.0
          d=5.0
          c=SQRT(1.0+d*bulkri)

          c_uwtr(i,j) = a/SQRT(1.0+2.0*b*bulkri/c)

        END IF

      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cucwtr

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CPTCWTR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cptcwtr(nx,ny,nz,ibgn,iend,jbgn,jend,zp,soiltyp,wspd,        &
           ptsfc,pt1,c_ptwtrneu,c_ptwtr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute c_ptwtr (the product of ustar and ptstar) at sea case.
!  Note:  a temperature scale defined as surface heat flux divided
!  by the friction velocity.
!
!  The quantity c_ptwtr is used by the subroutine SFCFLX to obtain
!  surface fluxes for both the unstable and stable cases.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong and X. Song
!  8/04/1994
!
!  MODIFICATION HISTORY:
!
!  2/27/95 (V.W. and X.S.)
!
!  1/12/96 (V.W. and X.S.)
!  Changed the calculation related to zo over the sea.
!  Added kvwtr to denote the Von Karman constant over the sea;
!  Set a lower limiter for zo, zolimit, and an upper limiter for z1, z1limit.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    iend     i-index where evaluation ends.
!    jend     j-index where evaluation ends.
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!    ptsfc    Potential temperature at the ground level (K)
!    pt1      Potential temperature at the 1st scalar point above
!             ground level, (K)
!
!
!  OUTPUT:
!
!    c_ptwtr  The product of ustar and ptstar. ptstar is temperature
!             scale (K), defined by surface heat flux / friction velocity
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

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins.
  INTEGER :: iend              ! i-index where evaluation ends.
  INTEGER :: jbgn              ! j-index where evaluation begins.
  INTEGER :: jend              ! j-index where evaluation ends.

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of staggered grid.
  INTEGER :: soiltyp(nx,ny)    ! Soil type at each point.

  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL :: ptsfc(nx,ny)         ! Potential temperature at the ground level
                               ! (K)
  REAL :: pt1   (nx,ny)        ! Potential temperature at the 1st scalar
                               ! point above ground level, (K)

  REAL :: c_ptwtrneu(nx,ny)       ! Product of ustar and ptstar
  REAL :: c_ptwtr(nx,ny)       ! Product of ustar and ptstar
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j

  REAL :: z1                   ! The height of 1st scalar point above
                               ! ground level (m)
  REAL :: zo                   ! Defined as surface momentum roughness
                               ! length
  REAL :: zt                   ! Defined as surface heat transfer
                               ! roughness length
  REAL :: dzo                  ! z1-zo
  REAL :: dzt                  ! z1-zt
  REAL :: z1dzo                ! z1/zo
  REAL :: z1dzt                ! z1/zt
  REAL :: xcdn                 ! Cdn (sea)
  REAL :: xcdh                 ! Hot-wired value of Cdh (sea)

  REAL :: bulkri               ! Bulk Richardson number
  REAL :: stabp                ! Monin-Obukhov STABility Parameter
                               ! (zeta)

  REAL :: y1,y0,psih           ! Intermediate parameters needed
  REAL :: qb3pb2
  REAL :: c7,c8
  REAL :: sb,qb,pb,thetab,tb   ! during  computations
  REAL :: a,b,c,d
  REAL :: tempan,sqrtqb
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  xcdh = 1.1E-3

  DO j=jbgn,jend
    DO i=ibgn,iend
      IF ( soiltyp(i,j) == 13) THEN

        xcdn = (0.4+0.07*wspd(i,j)) * 1.e-3
        z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
        z1  = MIN(z1,z1limit)
!
! calculate zo and zt
!
        zo =0.032*xcdn*wspd(i,j)*wspd(i,j)/9.8
        zo = MAX(zo,zolimit)
        zt = z1 * EXP( -kv*SQRT(xcdn)/(prantl0w*xcdh) )

        dzo = z1-zo
        z1dzo = z1/zo
        dzt = z1-zt
        z1dzt = z1/zt

        bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dzo*dzo/          &
                 (dzt*wspd(i,j)*wspd(i,j))

        IF (bulkri <= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  Unstable case: A modified formulation, which is similar to
!  equations (28)-(34) in Byun (1990).
!
!-----------------------------------------------------------------------
!
          bulkri = MAX (bulkri,-10.0)

          sb =bulkri/prantl0w

          qb=oned9*(c1w+c2w*sb*sb)
          pb=oned54*(c3w+c4w*sb*sb)

          qb3pb2=qb**3-pb*pb
          c7 = (z1*dzt/(dzo*dzo))*(ALOG(z1dzo)*ALOG(z1dzo)/ALOG(z1dzt))

          IF( qb3pb2 >= 0.0 ) THEN

            sqrtqb = SQRT(qb)
            tempan = MAX( -1.0, MIN( 1.0, pb/(sqrtqb**3) ) )

            thetab=ACOS(tempan)
            stabp =c7*(-2.0*SQRT(qb)*COS(thetab/3.0)+c5w)

          ELSE

            tb    =(SQRT(-qb3pb2)+ABS(pb))**oned3
            stabp =c7*(-(tb+qb/tb)+c5w)

          END IF
!
!-----------------------------------------------------------------------
!
!  According to a modified equation, which is similar to equation (14)
!  in Byun (1990).
!
!-----------------------------------------------------------------------
!
          c8=gammamw*stabp
          y1=SQRT(1.0 - c8)
          y0=SQRT(1.0 - c8/z1dzt)

          psih=2.0*ALOG((y1+1.0)/(y0+1.0))
!
!-----------------------------------------------------------------------
!
!  Compute ptstar via equation (11) in Byun (1981).
!
!-----------------------------------------------------------------------
!
          c_ptwtr(i,j)=kv / (prantl0w*(ALOG(z1dzt)-psih))

        ELSE
!
!-----------------------------------------------------------------------
!
!  Stable case: With the modified formulation in Louis et al (1981).
!
!-----------------------------------------------------------------------
!
          a=kv*kv/(prantl0w*LOG(z1dzo)*LOG(z1dzt))
          b=5.0
          d=5.0
          c=SQRT(1.0+d*bulkri)

          c_ptwtr(i,j) = SQRT(1.0+2.0*b*bulkri/c)
          c_ptwtr(i,j) = a*c_ptwtr(i,j)/(prantl0l*(1.0+3.0*b*bulkri*c))

        END IF
      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cptcwtr

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CUNEUWTR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cuneuwtr(nx,ny,nz,ibgn,iend,jbgn,jend,soiltyp,               &
           wspd, c_uneu)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute c_uneu (friction velocity) at the neutral state. The quantity
!  c_uneu is used by the subroutine SFCFLX to obtain surface fluxes.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong and X. Song
!  8/04/1994
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins.
!    iend     i-index where evaluation ends.
!    jbgn     j-index where evaluation begins.
!    jend     j-index where evaluation ends.
!
!    soiltyp  Soil type at each point
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!
!
!  OUTPUT:
!
!    c_uneu   Friction velocity at sea case.
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

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins.
  INTEGER :: iend              ! i-index where evaluation ends.
  INTEGER :: jbgn              ! j-index where evaluation begins.
  INTEGER :: jend              ! j-index where evaluation ends.

  INTEGER :: soiltyp(nx,ny)    ! Soil type at each point.
  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)


  REAL :: c_uneu   (nx,ny)     ! Frictional velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  DO j=jbgn,jend
    DO i=ibgn,iend

      IF (soiltyp(i,j) == 13) THEN
        c_uneu(i,j) = SQRT ((0.4+0.079*wspd(i,j)) * 1.e-3)
      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cuneuwtr
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE CPTNEUWTR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cptneuwtr(nx,ny,nz,ibgn,iend,jbgn,jend,                      &
           soiltyp,wspd,c_uneu,c_ptwtrneu)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong and X. Song
!  8/04/1994
!
!  MODIFICATION:
!
!  03/17/1997 (Yuhe Liu and Vince Wong)
!  Fixed a bug. Variable c_ptwtrneu was calculated in its inverse.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins
!    iend     i-index where evaluation ends
!    jbgn     j-index where evaluation begins
!    jend     j-index where evaluation ends
!
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!
!  OUTPUT:
!
!    c_ptwtrneu   Temperature scale (K), defined by
!             surface heat flux / friction velocity
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

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins.
  INTEGER :: iend              ! i-index where evaluation ends.
  INTEGER :: jbgn              ! j-index where evaluation begins.
  INTEGER :: jend              ! j-index where evaluation ends.

  INTEGER :: soiltyp (nx,ny)
  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL :: c_uneu (nx,ny)       !
  REAL :: c_ptwtrneu(nx,ny)    ! Temperature scale (K)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: xcdh
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  xcdh = 1.14E-3

  DO j=jbgn,jend
    DO i=ibgn,iend

      IF (soiltyp(i,j) == 13) THEN
        c_ptwtrneu(i,j) = xcdh/c_uneu(i,j)
      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cptneuwtr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE PBLDEPTH                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pbldepth(nx,ny,nz, u,v,w,ptprt,qv,ptbar,zp,j3, ptsflx,       &
           pbldpth)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calcualte time-dependent PBL depth.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong and X. Song
!  8/27/1994
!
!  Reference for option pbldopt=3 (not included yet):
!
!  1)" Simple Model of the Daytime Boundary Layer Height" by S-E.
!      Gryning and E. Batchvarova in 9th Symposium on Turbulence and
!      Diffusion, 379-382 (1990).
!  2) "A Rate Equation for the Nocturnal Boundary-Layer Height" by
!      F.T.M. Nieuwatadt and H. Tennekes in J. of Atmospheric Sciences
!      July 1981, pp. 1418-1428.
!
!  MODIFICATIONS:
!
!  3/6/1996 (M. Xue, V. Wong and Y. Liu)
!  Added option 2 for diagnostic calculatin of PBL depth.
!
!  06/06/96 (Ming Xue and Jinxing Zong)
!  PBL depth is now determined to be at the level where environmental
!  virtual potential temperate equals that at the first grid level.
!  Subroutine PBLDEPTH in sfcphy3d.f was modified.
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!
!    u        Total u-velocity (m/s)
!    v        Total v-velocity (m/s)
!    w        Total w-velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    qv       Specific humidity (kg/kg)
!
!    ptbar    Base state potential temperature (K)
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!    j3       Coordinate transformation Jacobian d(zp)/d(z)
!
!  OUTPUT:
!
!    pbldpth  Planetary boundary layer depth (m)
!
!  TEMPORARY:
!
!    ptv0     Work array for virtual pot. temp. at k=2 for scalar
!    ptv      Work array for virtual pot. temp.
!    parea    Work array for positive area
!    narea    Work array for negative area
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions
  INTEGER :: nt                ! Number of time levels
  PARAMETER( nt=3 )

  REAL :: u      (nx,ny,nz) ! Total u-velocity (m/s)
  REAL :: v      (nx,ny,nz) ! Total v-velocity (m/s)
  REAL :: w      (nx,ny,nz) ! Total w-velocity (m/s)
  REAL :: ptprt  (nx,ny,nz) ! Perturbation potential temperature (K)
  REAL :: qv     (nx,ny,nz) ! Specific humidity (kg/kg)

  REAL :: ptbar  (nx,ny,nz) ! Base state potential temperature (K)

  REAL :: zp     (nx,ny,nz) ! The physical height coordinate defined at
                            ! w-point of staggered grid.
  REAL :: j3     (nx,ny,nz) ! Coordinate transformation Jacobian  d(zp)/d(z)

  REAL :: ptsflx (nx,ny)    ! Surface heat flux (K*kg/(s*m**2))

  REAL :: pbldpth(nx,ny,nt) ! Planetary boundary layer depth (m)

  REAL :: ptv0   (nx,ny)    ! Temporary work array
  REAL :: ptv    (nx,ny)    ! Temporary work array
  REAL :: parea  (nx,ny)    ! Temporary work array
  REAL :: narea  (nx,ny)    ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL :: zptk0,ptvk0,zptk1,ptvk1,eps, amax,amin
  INTEGER :: pbldepth_set
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  IF ( pbldopt == 1 ) THEN

    DO j=1,ny-1
      DO i=1,nx-1
        pbldpth(i,j,2)=pbldpth0
      END DO
    END DO

  ELSE IF ( pbldopt == 2 ) THEN
!
!-----------------------------------------------------------------------
!
!  Diagnose the PBL depth based on the virtual potential temperature
!  profile. The PBL top is assumed to be at the level where the
!  environmental virtual potential temperature exceeds that at the
!  surface (first level above ground). If the atmosphere is stable
!  right above ground, the PBL depth is set to the thickness of the
!  layer below the first scalar point above ground.
!
!-----------------------------------------------------------------------
!
    DO j=1,ny-1
      DO i=1,nx-1
        ptv0(i,j)=(ptbar(i,j,2)+ptprt(i,j,2))*(1.+0.608*qv(i,j,2))
        parea(i,j) = 0.0
        narea(i,j) = 0.0
        pbldpth(i,j,2)=0.0  ! Set the initial depth to zero.
                            ! pbldpth also acts as a flag.
      END DO
    END DO

    eps = 1.0E-6
!
    DO k=3,nz-1

      DO j=1,ny-1
        DO i=1,nx-1

          IF( ABS(pbldpth(i,j,2)) < eps ) THEN
                              ! Check if pbldpth has been set.

            ptvk1=(ptbar(i,j,k  ) + ptprt(i,j,k  ))                     &
                       * (1. + 0.608*qv(i,j,k  ))
            ptvk0=(ptbar(i,j,k-1) + ptprt(i,j,k-1))                     &
                       * (1. + 0.608*qv(i,j,k-1))

            zptk0 = 0.5*(zp(i,j,k)+zp(i,j,k-1))
            zptk1 = 0.5*(zp(i,j,k+1)+zp(i,j,k))

            IF( ptvk1 > ptv0(i,j) ) THEN

              pbldpth(i,j,2) = zptk0 + (zptk1-zptk0)*                   &
                  (ptv0(i,j)-ptvk0)/(ptvk1-ptvk0) - zp(i,j,2)

            END IF

          END IF

        END DO
      END DO

    END DO

!-----------------------------------------------------------------------
!
!  When no stable layer is found above, the PBL top is at the model top.
!  The minimum PBL top height is at the first scalar point above ground.
!
!-----------------------------------------------------------------------

    pbldepth_set = 1

    DO j=1,ny-1
      DO i=1,nx-1
        IF(pbldpth(i,j,2) == 0.0) THEN
          pbldpth(i,j,2)=zp(i,j,nz-1)-zp(i,j,2)

          WRITE(6,'(/1x,a,a,2i4, 2(/1x,a,a))')                          &
              'Warning: The PBL is found to extend ',                   &
              'the entire model domain at i,j=',i,j,                    &
              'The atmosphere is either neutral or ',                   &
              'unstable throughout the domain depth.',                  &
              'In this case, PBL parameterization, if used, will operate', &
              ' on the entire domain depth.'
          pbldepth_set = 0
        END IF
        pbldpth(i,j,2) =                                                &
               MAX(pbldpth(i,j,2), 0.5*(zp(i,j,3)-zp(i,j,2)))
      END DO
    END DO

    IF(pbldepth_set == 0)THEN
      CALL a3dmax0(ptv0,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,amax,amin)
      WRITE(6,'(1x,2(a,e13.6))')'Max sfc virtual temperature =',amax
    END IF

  ELSE IF ( pbldopt == 3 ) THEN

    WRITE(6,'(/5x,a,i3,a,2(/5x,a))')                                    &
        'PBL depth calculation option=',pbldopt,' not avaiable.',       &
        'please reset parameter pbldopt in input file. ',               &
        'Program stopped in PBLDEPTH.'

    CALL arpsstop('arpsstop called from PBLDEPTH pbldopt incorrect ',1)

  END IF


  RETURN
END SUBROUTINE pbldepth

