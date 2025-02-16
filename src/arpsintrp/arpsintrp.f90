!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPSINTRP                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

PROGRAM arpsintrp
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This program interpolates gridded data from one ARPS grid to another.
!  It can be used to prepare data for running ARPS in a one-way nested
!  mode. It's exepcted to replace ARPSR2H in this capacity.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  3/27/1997. Written based on ARPSR2H and ARPSCVT.
!
!  MODIFICATION HISTORY:
!
!  7/31/1997 (M. Xue)
!  Added options for specifying the terrain for the new grid, in
!  addition to interpolated terrain.
!
!  4/24/1998 (D. Weber)
!  Added hydrostatic extrapolation calculations for pprt and pbar
!  when the fine grid terrain is lower in elevation than the
!  coarse grid terrain. (DO LOOP 1423)  Note: this update is
!  functional for bglopt=2 ONLY.
!
!  4/24/1998 (D. Weber)
!  Added option for quadratic interpolation in the vertical.
!      intrpvopt = 1 for linear
!      intrpvopt = 2 for quadratic
!
!  4/1/1999 (M. Xue)
!
!  Rewrote major portions of the program. Special treatment given
!  to grid points below input grid ground. When output grid
!  terrain falls below ground, base-state variables are reconstructed
!  from averages on horizontal planes to ensure x or y independency
!  below ground. This is done only for bglopt=4, however.
!  The new code also runs much faster.
!
!  2000/04/05 (Gene Bassett)
!  Added bglopt=5, which is similar to bglopt=4 but uses the
!  method in ext2arps for constructing the base sounding.
!
!  2000/04/17 (Ming Xue)
!  Added an option that allows one to specify input history data
!  at a constant time interval.
!
!  2000/05/20 (Gene Bassett)
!  Converted to F90, creating allocation and main subroutines.
!
!  2000/07/28 (Ming Xue)
!  Converted to F90 free format. Use ALLOCATABLE instead of
!  POINTER allocation to avoid double memory usage.
!
!  2000/10/31 (Gene Bassett)
!  Refined handling of interpolating soil variables.
!
!  2001/06/28 (Gene Bassett)
!  Added ntagopt, option to set surface winds in new grid to the surface
!  winds in the original grid for areas where new terrain is above the
!  original terrain.
!
!  1 June 2002 Eric Kemp
!  Soil variable updates.
!
!  13 June 2002 Eric Kemp
!  More updates, including improved processing of soil variables
!  from one grid to another.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'indtflg.inc'
  INCLUDE 'alloc.inc'
  INCLUDE 'exbc.inc'

  REAL, allocatable :: u     (:,:,:)  ! Total u-velocity (m/s).
  REAL, allocatable :: v     (:,:,:)  ! Total v-velocity (m/s).
  REAL, allocatable :: w     (:,:,:)  ! Total w-velocity (m/s).
  REAL, allocatable :: ptprt (:,:,:)  ! Perturbation potential temperature
                                      ! from that of base state atmosphere (Kelvin).
  REAL, allocatable :: pprt  (:,:,:)  ! Perturbation pressure from that
                                      ! of base state atmosphere (Pascal).
  REAL, allocatable :: qv    (:,:,:)  ! Water vapor specific humidity (kg/kg).
  REAL, allocatable :: qscalar(:,:,:,:)
  REAL, allocatable :: tke   (:,:,:)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, allocatable :: kmh   (:,:,:)  ! Horizontal turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, allocatable :: kmv   (:,:,:)  ! Vertical turb. mixing coef. for
                                      ! momentum. ( m**2/s )

  INTEGER, allocatable :: soiltyp(:,:,:)   ! Soil type
  REAL, allocatable ::    stypfrct(:,:,:)  ! Fraction of soil type
  INTEGER, allocatable :: vegtyp(:,:)      ! Vegetation type
  REAL, allocatable ::    lai    (:,:)     ! Leaf Area Index
  REAL, allocatable ::    roufns (:,:)     ! Surface roughness
  REAL, allocatable ::    veg    (:,:)     ! Vegetation fraction
  REAL, allocatable ::    ndvi   (:,:)     ! NDVI

  REAL, allocatable :: tsoil   (:,:,:,:)   ! soil temperature (K)
  REAL, allocatable :: qsoil   (:,:,:,:)   ! soil moisture
  REAL, allocatable :: wetcanp (:,:,:)     ! Canopy water amount
  REAL, allocatable :: snowdpth(:,:)       ! Snow depth (m)

  REAL, allocatable :: raing(:,:)          ! Grid supersaturation rain
  REAL, allocatable :: rainc(:,:)          ! Cumulus convective rain
  REAL, allocatable :: prcrate(:,:,:)      ! precipitation rate (kg/(m**2*s))
                                           ! prcrate(1,1,1) = total precip. rate
                                           ! prcrate(1,1,2) = grid scale precip. rate
                                           ! prcrate(1,1,3) = cumulus precip. rate
                                           ! prcrate(1,1,4) = microphysics precip. rate

  REAL, allocatable :: radfrc(:,:,:)       ! Radiation forcing (K/s)
  REAL, allocatable :: radsw (:,:)         ! Solar radiation reaching the surface
  REAL, allocatable :: rnflx (:,:)         ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet(:,:)       ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:)        ! Incoming longwave radiation

  REAL, allocatable :: usflx (:,:)         ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, allocatable :: vsflx (:,:)         ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, allocatable :: ptsflx(:,:)         ! Surface heat flux (K*kg/(m*s**2))
  REAL, allocatable :: qvsflx(:,:)         ! Surface moisture flux (kg/(m**2*s))

  REAL, allocatable :: ubar  (:,:,:)       ! Base state u-velocity (m/s).
  REAL, allocatable :: vbar  (:,:,:)       ! Base state v-velocity (m/s).
  REAL, allocatable :: wbar  (:,:,:)       ! Base state w-velocity (m/s).
  REAL, allocatable :: ptbar (:,:,:)       ! Base state potential temperature (K)
  REAL, allocatable :: pbar  (:,:,:)       ! Base state pressure (Pascal).
  REAL, allocatable :: rhobar(:,:,:)       ! Base state air density (kg/m**3)
  REAL, allocatable :: qvbar (:,:,:)       ! Base state water vapor specific humidity
                                           ! (kg/kg).
  REAL, allocatable :: x     (:)           ! The x-coord. of the physical and
                                           ! computational grid. Defined at u-point.
  REAL, allocatable :: y     (:)           ! The y-coord. of the physical and
                                           ! computational grid. Defined at v-point.
  REAL, allocatable :: z     (:)           ! The z-coord. of the computational grid.
                                           ! Defined at w-point on the staggered grid.
  REAL, allocatable :: zp    (:,:,:)       ! The physical height coordinate defined at
                                           ! w-point on the staggered grid.
  REAL, allocatable :: zpsoil(:,:,:)       ! The physical height coordinate defined at
                                           ! w-point on the staggered grid for soil model.
  REAL, allocatable :: hterain(:,:)        ! The height of terrain.

  REAL, allocatable ::    uprt   (:,:,:)   ! Perturbation u-velocity (m/s)
  REAL, allocatable ::    vprt   (:,:,:)   ! Perturbation v-velocity (m/s)
  REAL, allocatable ::    qvprt  (:,:,:)   ! Perturbation water vapor specific
                                           ! humidity (kg/kg)

  REAL, allocatable :: tem1  (:,:,:)       ! Temporary array
  REAL, allocatable :: tem2  (:,:,:)       ! Temporary array
  REAL, allocatable :: tem3  (:,:,:)       ! Temporary array
  REAL, allocatable :: tem4  (:,:,:)       ! Temporary array

! Temporary Soil arrays       WYH
  REAL, allocatable :: tsoiltem1(:,:,:), tsoiltem2(:,:,:)
  REAL, allocatable :: qsoiltem1(:,:,:), qsoiltem2(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Arrays on the new grid:
!
!-----------------------------------------------------------------------
!
  REAL, allocatable :: u1    (:,:,:)  ! Total u-velocity (m/s).
  REAL, allocatable :: v1    (:,:,:)  ! Total v-velocity (m/s).
  REAL, allocatable :: w1    (:,:,:)  ! Total w-velocity (m/s).
  REAL, allocatable :: ptprt1(:,:,:)  ! Perturbation potential temperature
                                      ! from that of base state atmosphere (Kelvin).
  REAL, allocatable :: pprt1 (:,:,:)  ! Perturbation pressure from that
                                      ! of base state atmosphere (Pascal).
  REAL, allocatable :: qv1   (:,:,:)  ! Water vapor specific humidity (kg/kg).
  REAL, allocatable :: qscalar1(:,:,:,:)
  REAL, allocatable :: tke1  (:,:,:)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, allocatable :: kmh1  (:,:,:)  ! Horizontal turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, allocatable :: kmv1  (:,:,:)  ! Vertical turb. mixing coef. for
                                      ! momentum. ( m**2/s )

  INTEGER, allocatable :: soiltyp1(:,:,:)   ! Soil type
  REAL, allocatable ::   stypfrct1(:,:,:)   ! Frction of soil type
  INTEGER, allocatable :: vegtyp1 (:,:)     ! Vegetation type
  REAL, allocatable ::    lai1    (:,:)     ! Leaf Area Index
  REAL, allocatable ::    roufns1 (:,:)     ! Surface roughness
  REAL, allocatable ::    veg1    (:,:)     ! Vegetation fraction

  REAL, allocatable :: tsoil1  (:,:,:,:)    ! soil temperature (K)
  REAL, allocatable :: qsoil1  (:,:,:,:)    ! soil moisture
  REAL, allocatable :: wetcanp1(:,:,:)      ! Canopy water amount
  REAL, allocatable :: snowdpth1(:,:)       ! Snow depth (m)

  REAL, allocatable :: raing1(:,:)          ! Grid supersaturation rain
  REAL, allocatable :: rainc1(:,:)          ! Cumulus convective rain
  REAL, allocatable :: prcrate1(:,:,:)      ! precipitation rate (kg/(m**2*s))
                                            ! prcrate(1,1,1) = total precip. rate
                                            ! prcrate(1,1,2) = grid scale precip. rate
                                            ! prcrate(1,1,3) = cumulus precip. rate
                                            ! prcrate(1,1,4) = microphysics precip. rate

  REAL, allocatable :: radfrc1(:,:,:)       ! Radiation forcing (K/s)
  REAL, allocatable :: radsw1 (:,:)         ! Solar radiation reaching the surface
  REAL, allocatable :: rnflx1 (:,:)         ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet1(:,:)       ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin1(:,:)        ! Incoming longwave radiation

  REAL, allocatable :: usflx1 (:,:)         ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, allocatable :: vsflx1 (:,:)         ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, allocatable :: ptsflx1(:,:)         ! Surface heat flux (K*kg/(m*s**2))
  REAL, allocatable :: qvsflx1(:,:)         ! Surface moisture flux (kg/(m**2*s))

  REAL, allocatable :: ubar1 (:,:,:)        ! Base state u-velocity (m/s).
  REAL, allocatable :: vbar1 (:,:,:)        ! Base state v-velocity (m/s).
  REAL, allocatable :: ptbar1(:,:,:)        ! Base state potential temperature (K)
  REAL, allocatable :: pbar1 (:,:,:)        ! Base state pressure (Pascal).
  REAL, allocatable :: rhobar1(:,:,:)       ! Base state air density (kg/m**3)
  REAL, allocatable :: qvbar1(:,:,:)        ! Base state water vapor specific humidity
                                            ! (kg/kg).

  REAL, allocatable :: x1    (:)            ! The x-coord. of the physical and
                                            ! computational grid. Defined at u-point.
  REAL, allocatable :: y1    (:)            ! The y-coord. of the physical and
                                            ! computational grid. Defined at v-point.
  REAL, allocatable :: z1    (:)            ! The z-coord. of the computational grid.
                                            ! Defined at w-point on the staggered grid.
  REAL, allocatable :: x1_out(:)
  REAL, allocatable :: y1_out(:)
  REAL, allocatable :: zp1   (:,:,:)        ! The physical height coordinate defined at
                                            ! w-point on the staggered grid.
  REAL, allocatable :: zpsoil1   (:,:,:)    ! The physical height coordinate defined at
                                            ! w-point on the staggered grid for soil model
  REAL, allocatable :: hterain1(:,:)        ! Terrain height (m)
  REAL, allocatable :: htrn1orig(:,:)        ! Terrain height (m)
  REAL, allocatable :: j11   (:,:,:)        ! Coordinate transformation Jacobian -d(zp)/dx.
  REAL, allocatable :: j21   (:,:,:)        ! Coordinate transformation Jacobian -d(zp)/dy.
  REAL, allocatable :: j31   (:,:,:)        ! Coordinate transformation Jacobian  d(zp)/dz.
  REAL, allocatable :: j3soil1(:,:,:)       ! Coordinate transformation Jacobian
                                            ! d(zpsoil)/dz.
  REAL, allocatable :: j3soilinv1(:,:,:)    ! Inverse of j3soil1

  REAL, allocatable :: tem11 (:,:,:)        ! Work array
  REAL, allocatable :: tem21 (:,:,:)        ! Work array

  REAL, allocatable :: wgtsx(:,:)  ! Weight for interpolation in x-dir for scalar points
  REAL, allocatable :: wgtsy(:,:)  ! Weight for interpolation in y-dir for scalar points
  REAL, allocatable :: wgtux(:,:)  ! Weight for interpolation in x-dir for u points
  REAL, allocatable :: wgtvy(:,:)  ! Weight for interpolation in y-dir for v points

  REAL, allocatable :: wgtz(:,:,:,:)  ! Weight for interpolation in z-dir for v points

  INTEGER, allocatable :: isx(:),jsy(:),iux(:),jvy(:),kz(:,:,:)

  REAL, allocatable :: zp1d1 (:)       ! Temporary array
  REAL, allocatable :: dzp1d1(:)       ! Temporary array

! Temporary Soil arrays       WYH
  REAL, allocatable :: tsoil1tem1(:,:,:), tsoil1tem2(:,:,:)
  REAL, allocatable :: qsoil1tem1(:,:,:), qsoil1tem2(:,:,:)
!
!-----------------------------------------------------------------------
!
!  More arrays.
!
!-----------------------------------------------------------------------
!
  REAL, allocatable :: xs (:), ys (:) ! x,y coord for scalar points
  REAL, allocatable :: xs1(:), ys1(:) ! x,y coord for scalar points

  REAL, allocatable :: xs_2d(:,:), ys_2d(:,:) ! used by subroutine extmnsnd

  REAL, allocatable :: temx1yz  (:,:,:)
  REAL, allocatable :: temx1y1z (:,:,:)
  REAL, allocatable :: temx1y1zb(:,:,:)

  REAL, allocatable :: temz1d1(:),temz1d2(:),temz1d3(:),temz1d4(:),     &
        temz1d5(:),temz1d6(:),temz1d7(:)

  REAL, allocatable :: zsnd(:),ptsnd(:),qvsnd(:),psnd(:),               &
        ubsnd(:),vbsnd(:)

  REAL, allocatable :: ptpsfc(:,:),ppsfc (:,:),qvsfc(:,:)
  REAL, allocatable :: qscalarsfc(:,:,:)
  REAL, allocatable :: tkesfc(:,:),kmhsfc(:,:),kmvsfc(:,:)
  REAL, allocatable :: ptbsfc(:,:),pbsfc(:,:),qvbsfc(:,:)
  REAL, allocatable :: htrnx1y1(:,:),radsfc(:,:)
  INTEGER, allocatable :: ktrnx1y1(:,:)

  INTEGER, allocatable :: ix(:,:),jy(:,:)
  REAL, allocatable :: xw(:,:),yw(:,:)

!
!-----------------------------------------------------------------------
!  nx, ny, nz: Dimensions of input grid.
!  nx1, ny1, nz1: Dimensions of output grid.
!-----------------------------------------------------------------------
!
  INTEGER :: nx       ! Number of input grid points in the x-direction
  INTEGER :: ny       ! Number of input grid points in the y-direction
  INTEGER :: nz       ! Number of input grid points in the z-direction
  INTEGER :: nzsoil   ! soil model levels

  INTEGER :: nx1      ! Number of output grid points in the x-direction
  INTEGER :: ny1      ! Number of output grid points in the y-direction
  INTEGER :: nz1      ! Number of output grid points in the z-direction
  INTEGER :: nzsoil1  ! soil model levels

  INTEGER :: nxyz,nxy,nxyz1,nxy1
  INTEGER :: nxyzsoil, nxyzsoil1
!
!-----------------------------------------------------------------------
!
!  Soil types.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nstypin           ! Number of soil types in input data
  INTEGER :: nstyp1            ! Number of soil types for ouput data
!
!-----------------------------------------------------------------------
!
!  Misc
!
!-----------------------------------------------------------------------
!
  INTEGER :: lvlprof           ! Number of levels in 1-d average sounding
!
!-----------------------------------------------------------------------
!
  INTEGER :: hinfmt,nhisfile,lengbf,nf,lenfil
  INTEGER, PARAMETER  :: nhisfile_max = 200
  CHARACTER (LEN=256) :: grdbasfn, hisfile(nhisfile_max)
  INTEGER :: ireturn

  NAMELIST /output_dims/ nx1,ny1,nz1

  REAL :: dzsoil1,zrefsoil1

! Commented out variables already declared in grid.inc.  But keep in
! the namelist!!!
!
!  INTEGER :: soilstrhopt
!  REAL :: soildzmin,soildlayer1,soildlayer2,soilstrhtune

  NAMELIST /newgrid_soil/ soilmodel_option,nzsoil1,dzsoil1,zrefsoil1, &
                          soilstrhopt,soildzmin,soildlayer1, &
                          soildlayer2,soilstrhtune

  INTEGER, PARAMETER :: max_plevels=1000

  REAL :: plevels(max_plevels)
                       ! Array contain the values of pressure levels to
                       ! which fields will be interpolated.
!
!-----------------------------------------------------------------------
!
!  Parameters for the output grid.
!
!  Note:
!
!  Given nx1, ny1 and nz1, the physical domain size of the refined
!  grid will be xl1=(nx1-3)*dx1 by yl1=(ny1-3)*dy1 by zh1=(nz1-3)*dz1.
!  Dx1, dy1 and dz1 are the grid intervals of the refined grid.
!
!-----------------------------------------------------------------------
!
  REAL :: xorig1, yorig1, zorig1
  REAL :: xctr1 , yctr1
  REAL :: ctrlat1, ctrlon1
  REAL :: origlat, origlon

  REAL :: dx1,dy1     ! Grid intervals of the refined grid.

  INTEGER :: strhopt1 ! Vertical grid stretching option.
                      ! = 0, no stretching in vertical.
                      ! >= 1, with stretching in vertical.
  REAL :: dz1         ! Average grid spacing in vertical direction in
                      ! transformed computational space (m).
  REAL :: dzmin1      ! Minimun grid spacing in vertical direction in
                      ! physcal space (m).

  REAL :: zrefsfc1    ! The reference height of the surface
                      ! (ground level) (m)

  REAL :: dlayer11    ! The depth of the lower layer with uniform
                      ! (dz=dzmin) vertical spacing (m)

  REAL :: dlayer21    ! The depth of the mid layer with stetched
                      ! vertical spacing (m)

  REAL :: strhtune1   ! Tuning parameter for stretching option 2
                      ! A Value between 0.2 and 5.0 is appropriate.
                      ! A larger value gives a more linear stretching.

  REAL :: zflat1      ! The height at which the grid levels
                      ! becomes flat in the terrain-following
                      ! coordinate transformation (m).
!
!-----------------------------------------------------------------------
!
!  Terrain option parameters for the new grid:
!
!-----------------------------------------------------------------------
!

  INTEGER :: ternopt1        ! Model terrain option.
                             ! = 0, no terrain, the ground is flat;
                             ! = 1, bell-shaped mountain;
                             ! = 2, terrain data read from terrain data
                             !      base (not implemented yet)
  INTEGER :: ntmerge         ! Number of zones over which to merge original
                             ! and new terrain at the bondaries
                             ! (used when ternopt1=4).
  INTEGER :: ternfmt1        ! Terrain data file format.
  INTEGER :: mntopt1         ! Option for choosing idealized mountain
                             ! type.
  REAL :: hmount1         ! The mountain height (m)
  REAL :: mntwidx1        ! The half-width of bell-shaped mountain
                          ! in x-dir.
  REAL :: mntwidy1        ! The half-width of bell-shaped mountain
                          ! in y-dir.
  REAL :: mntctrx1        ! The x coordinate of the bell-shaped
                          ! mountain center.
  REAL :: mntctry1        ! The y coordinate of the bell-shaped
                          ! mountain center.
  CHARACTER (LEN=256) :: terndta1  ! Name of the terrain data file

  REAL :: zflat11,za,zb
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nunit
  INTEGER :: i,j,k,nq
  REAL    :: amin, amax
  REAL    :: tmp1,tmp2
  REAL    :: radnd,pi2,z0
  INTEGER :: k0,kk
  CHARACTER (LEN=256) :: basdmpfn
  INTEGER             :: lbasdmpf
  CHARACTER (LEN=256) :: ternfn,sfcoutfl,soiloutfl,temchar
  INTEGER             :: lternfn,lfn
  INTEGER :: iss

  DOUBLE PRECISION :: ntmergeinv, mfac
  INTEGER :: idist
  INTEGER :: houtfmt, nchin, nchout
  CHARACTER (LEN=256) :: filename

  INTEGER :: grdbas

  REAL :: time
  INTEGER :: vroutcnt
  DATA    vroutcnt /0/

  INTEGER :: nfile, length, lenstr
  CHARACTER (LEN=132) :: timsnd
  CHARACTER (LEN=80) :: new_runname
  INTEGER :: tmstrln
  REAL :: xeps, yeps
  REAL :: ctrx,ctry,swx,swy,alatpro(2),sclf,dxscl,dyscl
  INTEGER :: bglopt,ntagopt
  REAL :: misvalue,aghght

  REAL :: snddelz,delz,dtdz,tbark1,tbark,ttotk1,ttotk,lnpbar,lnptot

  INTEGER :: i1mn,i1mx,j1mn,j1mx

  INTEGER :: npoint_below_ground ! flag indicating if any point in the
                                 ! output grid terrain is below input grid terrain
  INTEGER :: redo_base_state
  REAL :: zpmin,zpmax,zpmin1,pbartop, tem,tema,a,b,c,d

  INTEGER :: realtime,ntries,sleeptime
!
!-----------------------------------------------------------------------
!
!  namelist Declarations:
!
!-----------------------------------------------------------------------
!
  NAMELIST /jobname/ runname

  INTEGER :: z_or_p   ! Control parameter denoting if the output grid is
                      ! in height or pressure coordinates.
  INTEGER :: xy_or_ll ! Control parameter denoting if (xctr1,yctr1) or
                      ! (ctrlat1,ctrlon1) are to be used to specify the new
                      ! grid center.
  INTEGER :: intrphopt  ! Option for horizontal interpolation
                        ! = 1 linear, =2 quadratic
  INTEGER :: intrpvopt  ! Option for vertical interpolation
                        ! = 1 linear, =2 quadratic

  INTEGER :: snap_to_grid, same_res
  INTEGER :: istatus, ib,jb

  INTEGER, PARAMETER :: MAX_GRD=200
  INTEGER :: noutgrds
  REAL :: xctr_grd(MAX_GRD), yctr_grd(MAX_GRD)
  REAL :: clat_grd(MAX_GRD), clon_grd(MAX_GRD)
  CHARACTER (LEN=80 ) :: name_grd(MAX_GRD)
  INTEGER :: ng

  REAL ::  ctrlat_sv, ctrlon_sv, latitud_sv, dx_sv, dy_sv, dz_sv
  REAL :: dzmin_sv, zrefsfc_sv, dlayer1_sv, dlayer2_sv, zflat_sv, strhtune_sv
  INTEGER :: strhopt_sv, nstyp_sv

  NAMELIST /newgrid/ z_or_p,xy_or_ll,strhopt1,xctr1,yctr1,              &
          ctrlat1,ctrlon1,snap_to_grid,same_res,dx1,dy1,dz1,dzmin1,     &
          zrefsfc1,dlayer11,dlayer21,strhtune1,zflat1,plevels,nstyp1,   &
          noutgrds,xctr_grd,yctr_grd,clat_grd,clon_grd,name_grd

  NAMELIST /newterrain/ ternopt1,mntopt1,hmount1,mntwidx1,mntwidy1,     &
            mntctrx1,mntctry1,terndta1,ternfmt1,ntmerge

  NAMELIST /bgloption/ bglopt,misvalue,intrphopt,intrpvopt,ntagopt,aghght

  NAMELIST /sfc_data/ sfcdat,styp,vtyp,lai0,roufns0,veg0,               &
            sfcdtfl,sfcfmt

  NAMELIST /output/ dirname,exbcdmp,exbchdfcompr,hdmpfmt,grbpkbit,hdfcompr,&
            grdout,basout,varout,mstout,rainout,prcout,iceout,          &
            tkeout, trbout,sfcout,landout,radout,flxout,                &
            qcexout,qrexout,qiexout,qsexout,qhexout,qgexout,nqexout,zqexout, &
            totout,filcmprs,readyfl,sfcdmp,soildmp,terndmp,ngbrz,zbrdmp
  NAMELIST /process/ realtime,ntries,sleeptime
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: cpu0, gcpu_time, f_cputime
  REAL :: xctr1_old,yctr1_old,xsw,ysw,xsw1,ysw1,xsw1_old,ysw1_old
  REAL :: ctrlat1_old, ctrlon1_old

  REAL :: rhmax = 1.0
  REAL :: qvsat

  REAL :: hght0, dhght, fac
  REAL :: ptop

  INTEGER :: nstyps_sfcdat
  INTEGER, ALLOCATABLE :: soiltyp_his(:,:,:)   ! Soil type
  REAL, ALLOCATABLE :: wetcanp_his(:,:,:)      ! Wet canopy
  REAL, ALLOCATABLE :: tsoil_his(:,:,:,:)      ! Soil temperature
  REAL, ALLOCATABLE :: qsoil_his(:,:,:,:)      ! Soil moisture
  REAL, ALLOCATABLE :: x2d1(:,:),y2d1(:,:)
  INTEGER, ALLOCATABLE :: i2d1(:,:),j2d1(:,:),k3d1(:,:,:)
  INTEGER :: ii,is
  REAL :: dzsoil_tmp
  REAL :: frctot

  REAL, allocatable :: dxfld(:)
  REAL, allocatable :: dyfld(:)
  REAL, allocatable :: rdxfld(:)
  REAL, allocatable :: rdyfld(:)
  REAL, ALLOCATABLE :: rdzsoilfld(:,:,:)
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
  cpu0 = f_cputime()

  WRITE(6,'(/9(/2x,a)/)')                                               &
     '###############################################################', &
     '###############################################################', &
     '###                                                         ###', &
     '###                Welcome to ARPSINTRP                     ###', &
     '###      This program converts the history dump data        ###', &
     '###      sets generated by ARPS, between various formats.   ###', &
     '###                                                         ###', &
     '###############################################################', &
     '###############################################################'

!
!-----------------------------------------------------------------------
!
!  Read in the input parameters.
!
!-----------------------------------------------------------------------
!
  CALL get_input_file_names(5,hinfmt,grdbasfn,hisfile,nhisfile)

  lenfil = len_trim(hisfile(1))

  CALL get_dims_from_data(hinfmt,hisfile(1)(1:lenfil),                    &
       nx,ny,nz,nzsoil,nstypin, ireturn)

  nstypin = max(1,nstypin)

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(4(a,i5))')  &
          'nx =',nx,', ny=',ny,', nz=',nz,', nzsoil',nzsoil

  READ(5,output_dims, END=105)
  WRITE(6,'(/a,a/)')                                                    &
      'NAMELIST block output_dims successfully read.'

  READ(5,newgrid_soil, END=105)
  WRITE(6,'(/a,a/)')                                                    &
      'NAMELIST block newgrid_soil successfully read.'

  WRITE(6,'(4(a,i5))')  &
          'nx1=',nx1,', ny1=',ny1,', nz1=',nz1,', nzsoil1',nzsoil1

  GO TO 10

  105   WRITE(6,'(/a,a/)')                                              &
          'Error reading NAMELIST block intrp_dims. ',                  &
          'Program ARPSINTRP stopped.'
  STOP 1

  10    CONTINUE

  WRITE (6,'(/a,g13.3/)') "Current memory allocation (in words):",      &
               current_memory_use


  IF (max_plevels < nz1) THEN
    WRITE (6,*) "ARPSINTRP: ERROR, nz1 < max_plevels.  Stopping."
    WRITE (6,*) "nz1 =",nz1,"  max_plevels = ",max_plevels
    STOP 1
  END IF

  mgrid = 1
  nestgrd = 0
!
!-----------------------------------------------------------------------
!
!  Set the default parameters
!
!-----------------------------------------------------------------------
!
  runname   = 'intrp'

  z_or_p    = 1
  xy_or_ll  = 1
  xctr1     = 512000.0
  yctr1     = 512000.0
  ctrlat1   = 34.0
  ctrlon1   = -98.0
  dx1       =  4000.0
  dy1       =  4000.0
  dz1       =  350.0
  strhopt1  = 2
  dzmin1    = 20.0
  zrefsfc1  = 0.0
  dlayer11  = 0.0
  dlayer21  = 1.0E5
  strhtune1 = 1.0
  zflat1    = 1.0E5
  nstyp1    = -1

  noutgrds = 0
  xctr_grd = 0
  yctr_grd = 0
  clat_grd = ctrlat1
  clon_grd = ctrlon1
  name_grd = 'NULL'

  ternopt1  = 3
  ntmerge = 1
  mntopt1   = 1
  hmount1   =     0.000
  mntwidx1  = 10000.000
  mntwidy1  = 10000.000
  mntctrx1  = 10000.000
  mntctry1  = 10000.000
  terndta1  = 'arpstern.dat'
  ternfmt1  = 1

  bglopt    = 1
  misvalue  = -9999.0
  intrphopt = 1
  intrpvopt = 1
  ntagopt = 0
  aghght = 500.0

  sfcdat  = 1
  styp    = 3
  vtyp    = 10
  lai0    = 0.31
  roufns0 = 0.1
  veg0    = 0.0
  sfcdtfl = 'arpssfc.data'
  sfcfmt  = 1

  dirname   = './'
  exbcdmp   = 1
  exbchdfcompr = 0
  hdmpfmt   = 1
  grbpkbit  = 16
  hdfcompr  = 0
  filcmprs  = 0
  readyfl   = 1
  basout    = 0
  grdout    = 0
  varout    = 1
  mstout    = 1
  iceout    = 1
  tkeout    = 1
  trbout    = 0
  rainout   = 0
  sfcout    = 0
  snowout   = 0
  landout   = 0
  qcexout   = 0
  qrexout   = 0
  qiexout   = 0
  qsexout   = 0
  qhexout   = 0
  qgexout   = 0
  nqexout   = 0
  zqexout   = 0

  ngbrz = 5
  zbrdmp = 10000.0

  sfcdmp    = 1
  soildmp   = 1
  terndmp   = 1
  same_res  = 0
  snap_to_grid = 0

  realtime  = 0
  ntries    = 240
  sleeptime = 30

  READ (5,jobname,END=100)
  WRITE(6,'(/a/)') 'Sucessfully read namelist block JOBNAME.'

  WRITE(6,'(/2x,a,a)') 'The name of this run is: ', runname
  new_runname = runname
!
!-----------------------------------------------------------------------
!
!  Set the output grid and the variable control parameters
!
!-----------------------------------------------------------------------
!
  READ (5,newgrid)
  WRITE(6,'(/a/)') 'Sucessfully read namelist block NEWGRID.'

  PRINT*
  PRINT*,' Input parameters for the new refined grid:'
  PRINT*

  PRINT*,' z_or_p   = ', z_or_p
  PRINT*,' xy_or_ll = ', xy_or_ll
  PRINT*,' dx1      = ',dx1
  PRINT*,' dy1      =  ',dy1
  PRINT*,' strhopt1 = ',strhopt1
  PRINT*,' dz1      =  ',dz1
  PRINT*,' dzmin1   =  ',dzmin1
  PRINT*,' xctr1    =  ',xctr1
  PRINT*,' yctr1    =  ',yctr1
  PRINT*,' ctrlat1  = ',ctrlat1
  PRINT*,' ctrlon1  = ',ctrlon1
  PRINT*,' plevels  = ', (plevels(k),k=1,nz1)
  PRINT*,' nstyp1   = ',nstyp1

  IF (z_or_p == 2) THEN
    k = 1
    DO WHILE (.TRUE.)
      IF (plevels(k) < 0.01) THEN    ! plevels must > 0
        k = k-1
        EXIT
      END IF
      IF (plevels(k+1) > plevels(k)) THEN
        WRITE(6,'(/1x,a/)') 'ERROR: variable "plevels" must be in decending order.'
        STOP
      END IF
      k = k + 1
    END DO

    IF (k /= nz1) THEN
      nz1 = k
      WRITE(6,'(/,1x,a,I2,a,I2,a,/)') 'WARNING: nz1 is reset to ',k,    &
         ', since only ',k,' valid values are found in variable plevels.'
      WRITE(6,*) ' plevels = ',(plevels(k),k=1,nz1)
      WRITE(6,*)
    END IF

    IF (nz1 < 1) THEN
      WRITE(6,'(1x,a,I3,a)') 'ERROR: Wrong size of nz1 = ',nz1,'.'
      STOP
    END IF
  END IF

  PRINT*,'noutgrds = ',noutgrds
  IF (noutgrds > MAX_GRD) THEN
    WRITE (*,*) "ERROR, noutgrds greater than MAX_GRD. ",  &
       "Set noutgrds to less than or equal to ",MAX_GRD," ABORTING!"
    STOP
  ENDIF
  IF (noutgrds > 0) THEN
    DO ng = 1,noutgrds
      PRINT*,'xctr_grd(',ng,') = ',xctr_grd(ng)
      PRINT*,'yctr_grd(',ng,') = ',yctr_grd(ng)
      PRINT*,'clat_grd(',ng,') = ',clat_grd(ng)
      PRINT*,'clon_grd(',ng,') = ',clon_grd(ng)
      PRINT*,'name_grd(',ng,') = ',name_grd(ng)
    END DO
  ENDIF

!
!-----------------------------------------------------------------------
!
!  Specify the terrain initialization option parameters:
!
!-----------------------------------------------------------------------
!

  READ (5,newterrain,END=100)
  WRITE(6,'(/a)') 'Sucessfully read namelist block NEWTERRAIN.'

  WRITE(6,'(1x,a,i4)') 'ternopt1 =', ternopt1
  WRITE(6,'(1x,a,i4)') 'ntmerge  =', ntmerge
  WRITE(6,'(1x,a,i4)') 'mntopt1  =', mntopt1
  WRITE(6,'(1x,a,f10.3)') 'hmount1 =',hmount1
  WRITE(6,'(1x,a,f10.3)') 'mntwidx1=',mntwidx
  WRITE(6,'(1x,a,f10.3)') 'mntwidy1=',mntwidy
  WRITE(6,'(1x,a,f10.3)') 'mntctrx1=',mntctrx
  WRITE(6,'(1x,a,f10.3)') 'mntctry1=',mntctry

  lenstr = 256
  CALL strlnth( terndta1, lenstr)
  WRITE(6,'(2x,a,a/)') 'terndta1 = ',terndta1(1:lenstr)
  WRITE(6,'(1x,a,i4)') 'ternfmt1 =', ternfmt1
  ternfmt = ternfmt1


  READ(5,bgloption,END=100)
  WRITE(6,'(/a)') 'Sucessfully read namelist block BLGOPTION.'

  WRITE(6,'(1x,a,i4)')    'bglopt  =', bglopt
  WRITE(6,'(1x,a,f10.3)') 'misvalue=', misvalue
  WRITE(6,'(1x,a,i4)')    'intrphopt  =', intrphopt
  WRITE(6,'(1x,a,i4)')    'intrpvopt  =', intrpvopt
  WRITE(6,'(1x,a,i4)')    'ntagopt  =', ntagopt
  WRITE(6,'(1x,a,f10.3)') 'aghght  =', aghght

!
!-----------------------------------------------------------------------
!
!  Set the control parameters for soil data:
!
!-----------------------------------------------------------------------
!
  READ (5,sfc_data,END=100)
  WRITE(6,'(/a/)') 'Successfully read namelist block SOIL_DATA.'

  IF (sfcdat == 1) THEN ! Sfc characteristics set in input file.
    nstyp = 1
    nstypin = nstyp
  ELSE IF (sfcdat == 2 .OR. sfcdat == 3) THEN
    ! Soiltype(s) etc. set in sfcdat file.

    CALL get_nstyps_from_sfcdat(nstyps_sfcdat,sfcdtfl)
    nstyp = MAX(nstypin,nstyps_sfcdat)
  ELSE ! Use soiltype etc. in history dump.
    nstyp = nstypin
  ENDIF

  IF (nstyp1 .eq. -1) nstyp1 = nstyp

  WRITE (6,'(a,i3)') "Number of soil types to be used for input data:",nstypin
!
!-----------------------------------------------------------------------
!
!  Set the control parameters for output:
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(a)')                                                        &
      ' Reading in control parameters for the output data files..'

  READ (5,output,END=100)
  WRITE(6,'(a)') 'Successfully read namelist block OUTPUT.'

  houtfmt = hdmpfmt

!  IF( houtfmt.eq.0 ) THEN
!    write(6,'(/1x,a/)') 'Output format is 0, no data is dumped.'
!    STOP 9103
!  ENDIF

  IF ( houtfmt == 10 .AND. grbpkbit <= 0 ) THEN
    WRITE(6,'(a,a,i2/a)')                                               &
        'The bit width for packing GRIB data was invalid, ',            &
        'The old value was ', grbpkbit, 'Reset it to 16 bits'
    grbpkbit = 16
  END IF

  totout = 1

  IF (sfcout == 1) THEN
    snowout = 1
  END IF

!
!-----------------------------------------------------------------------
!
!  Set the control parameters for process:
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(a)')                                                        &
      ' Reading in control parameters for the process data files..'

  READ (5,process,END=100)
  WRITE(6,'(a)') 'Successfully read namelist block PROCESS.'

  IF ( realtime /= 0 ) THEN
	write(6,*) 'Realtime mode'
  ENDIF

  GO TO 15

  100   WRITE(6,'(a)')                                                  &
          'Error reading NAMELIST file. Program ARPSINTRP stopped.'
  STOP 9104

  15    CONTINUE

  ldirnam=LEN(dirname)
  CALL strlnth( dirname , ldirnam)

  lengbf=len_trim(grdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', grdbasfn(1:lengbf)
!
!-----------------------------------------------------------------------
!
!  Allocate arrays.
!
!-----------------------------------------------------------------------
!
  WRITE (6,'(/a,g13.3/)') "Current memory allocation (in words):",      &
               current_memory_use

!
  lvlprof=MAX(601,nz)
  nxyz = nx*ny*nz
  nxyzsoil = nx*ny*nzsoil
  nxy = nx*ny

  allocate(u(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'u')
  allocate(v(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'v')
  allocate(w(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'w')
  allocate(ptprt(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'ptprt')
  allocate(pprt(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'pprt')

  allocate(qv(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'qv')
  allocate(qscalar(nx,ny,nz,nscalar))
  CALL alloc_status_accounting(istatus,nxyz*nscalar,'qscalar')
  allocate(tke(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'tke')
  allocate(kmh(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'kmh')
  allocate(kmv(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'kmh')

  ALLOCATE(soiltyp_his(nx,ny,nstyp),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy*nstyp,'soiltyp_his')
  ALLOCATE(wetcanp_his(nx,ny,0:nstyp),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy*(nstyp+1),'wetcanp_his')
  ALLOCATE(tsoil_his(nx,ny,nzsoil,0:nstyp),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy*nzsoil*(nstyp+1),'tsoil_his')
  ALLOCATE(qsoil_his(nx,ny,nzsoil,0:nstyp),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy*nzsoil*(nstyp+1),'qsoil_his')

  allocate(soiltyp(nx,ny,nstyp),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy*nstyp,'soiltyp')
  allocate(stypfrct(nx,ny,nstyp),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy*nstyp,'stypfrct')
  allocate(vegtyp(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'vegtyp ')
  allocate(lai(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'lai    ')
  allocate(roufns(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'roufns ')
  allocate(veg(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'veg    ')
  allocate(ndvi(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'ndvi   ')

  allocate(tsoil(nx,ny,nzsoil,0:nstyp),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy*nzsoil*(nstyp+1),'tsoil  ')
  allocate(qsoil(nx,ny,nzsoil,0:nstyp),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy*nzsoil*(nstyp+1),'qsoil  ')
  allocate(wetcanp(nx,ny,0:nstyp),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy*(nstyp+1),'wetcanp')
  allocate(snowdpth(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'snowdpth')

  allocate(raing(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'raing')
  allocate(rainc(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'rainc')
  allocate(prcrate(nx,ny,4),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy*4,'prcrate')

  allocate(radfrc(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'radfrc')
  allocate(radsw(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'radsw ')
  allocate(rnflx(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'rnflx ')
  allocate(radswnet(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'radswnet')
  allocate(radlwin(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'radlwin')
  allocate(usflx(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'usflx ')
  allocate(vsflx(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'vsflx ')
  allocate(ptsflx(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'ptsflx')
  allocate(qvsflx(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'qvsflx')

  allocate(ubar(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'ubar')
  allocate(vbar(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'vbar')
  allocate(wbar(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'wbar')
  allocate(ptbar(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'ptbar')
  allocate(pbar(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'pbar')
  allocate(rhobar(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'rhobar')
  allocate(qvbar(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'qvbar')

  allocate(x(nx),stat=istatus)
  CALL alloc_status_accounting(istatus,nx,'x')
  allocate(y(ny),stat=istatus)
  CALL alloc_status_accounting(istatus,ny,'y')
  allocate(z(nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nz,'z')
  allocate(zp(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'zp')
  allocate(zpsoil(nx,ny,nzsoil),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyzsoil,'zpsoil')
  allocate(hterain(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'hterain')

  allocate(uprt(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'uprt')
  allocate(vprt(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'vprt')
  allocate(qvprt(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'qvprt')

  allocate(tem1(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'tem1')
  allocate(tem2(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'tem2')
  allocate(tem3(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'tem3')
  allocate(tem4(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'tem4')

  allocate(u1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'u1')
  allocate(v1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'v1')
  allocate(w1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'w1')
  allocate(ptprt1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'ptprt1')
  allocate(pprt1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz,'pprt1')

  nxyz1 = nx1*ny1*nz1
  nxyzsoil1 = nx1*ny1*nzsoil1
  nxy1  = nx1*ny1

  print*,'Start allocating arrays for output grid.'

  ALLOCATE(x2d1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*ny1,'x2d1')
  ALLOCATE(y2d1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*ny1,'y2d1')
  ALLOCATE(i2d1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*ny1,'i2d1')
  ALLOCATE(j2d1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*ny1,'j2d1')
  ALLOCATE(k3d1(nx1,ny1,nzsoil1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*ny1*nzsoil1,'k3d1')
  ALLOCATE(dxfld(nx),stat=istatus)
  CALL alloc_status_accounting(istatus,nx,'dxfld')
  dxfld = 0.0
  ALLOCATE(rdxfld(nx),stat=istatus)
  CALL alloc_status_accounting(istatus,nx,'rdxfld')
  rdxfld = 0.0
  ALLOCATE(dyfld(ny),stat=istatus)
  CALL alloc_status_accounting(istatus,ny,'dyfld')
  dyfld = 0.0
  ALLOCATE(rdyfld(ny),stat=istatus)
  CALL alloc_status_accounting(istatus,ny,'rdyfld')
  rdyfld = 0.0
  ALLOCATE(rdzsoilfld(nx,ny,nzsoil),stat=istatus)
  CALL alloc_status_accounting(istatus,nx*ny*nzsoil,'rdzsoilfld')
  rdzsoilfld = 0.0

  allocate(qv1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'qv1')
  qv1 = 0.0
  allocate(qscalar1(nx1,ny1,nz1,nscalar),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1*nscalar,'qscalar1')
  qscalar1 = 0.0
  allocate(tke1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'tke1')
  tke1= 0.0
  allocate(kmh1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'kmh1')
  kmh1= 0.0
  allocate(kmv1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'kmv1')
  kmv1= 0.0

  allocate(soiltyp1(nx1,ny1,nstyp1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1*nstyp1,'soiltyp1')
  allocate(stypfrct1(nx1,ny1,nstyp1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1*nstyp1,'stypfrct1')
  soiltyp1 = 0
  stypfrct1 = 0.0
  allocate(vegtyp1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'vegtyp1')
  vegtyp1 = 0
  allocate(lai1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'lai1')
  lai1 = 0.0
  allocate(roufns1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'roufns1')
  roufns1 = 0.0
  allocate(veg1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'veg1')
  veg1 = 0.0

  allocate(tsoil1(nx1,ny1,nzsoil1,0:nstyp1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1*nzsoil1*(nstyp1+1),'tsoil1')
  tsoil1 = 0.0
  allocate(qsoil1(nx1,ny1,nzsoil1,0:nstyp1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1*nzsoil1*(nstyp1+1),'qsoil1')
  qsoil1 = 0.0
  allocate(wetcanp1(nx1,ny1,0:nstyp1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1*(nstyp1+1),'wetcanp1')
  wetcanp1 = 0.0
  allocate(snowdpth1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'snowdpth1')
  snowdpth1 = 0.0

  allocate(raing1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'raing1')
  raing1 = 0.0
  allocate(rainc1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'rainc1')
  rainc1 = 0.0
  allocate(prcrate1(nx1,ny1,4),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1*4,'prcrate1')
  prcrate1 = 0.0

  allocate(radfrc1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'radfrc1')
  radfrc1 = 0.0
  allocate(radsw1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'radsw1')
  radsw1 = 0.0
  allocate(rnflx1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'rnflx1')
  rnflx1 = 0.0
  allocate(radswnet1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'radswnet1')
  radswnet1 = 0.0
  allocate(radlwin1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'radlwin1')
  radlwin1 = 0.0

  allocate(usflx1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'usflx1')
  usflx1 = 0.0
  allocate(vsflx1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'vsflx1')
  vsflx1 = 0.0
  allocate(ptsflx1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'ptsflx1')
  ptsflx1 = 0.0
  allocate(qvsflx1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'qvsflx1')
  qvsflx1 = 0.0

  allocate(ubar1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'ubar1')
  ubar1 = 0.0
  allocate(vbar1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'vbar1')
  vbar1 = 0.0
  allocate(ptbar1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'pbar1')
  ptbar1 = 0.0
  allocate(pbar1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'ptbar1')
  pbar1 = 0.0
  allocate(rhobar1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'rhobar1')
  rhobar1 = 0.0
  allocate(qvbar1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'qvbar1')
  qvbar1 = 0.0

  allocate(x1(nx1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1,'x1')
  x1 = 0.0
  allocate(y1(ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,ny1,'y1')
  y1 = 0.0
  allocate(x1_out(nx1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1,'x1_out')
  x1_out = 0.0
  allocate(y1_out(ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,ny1,'y1_out')
  y1_out = 0.0
  allocate(z1(nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nz1,'z1')
  z1 = 0.0
  allocate(zp1(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'zp1')
  zp1 = 0.0
  allocate(zpsoil1(nx1,ny1,nzsoil1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyzsoil1,'zpsoil1')
  zpsoil1 = 0.0
  allocate(hterain1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'hterain1')
  hterain1 = 0.0
  allocate(htrn1orig(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'htrn1orig')
  htrn1orig = 0.0
  allocate(j11(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'j11')
  j11 = 0.0
  allocate(j21(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'j21')
  j21 = 0.0
  allocate(j31(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'j31')
  j31 = 1.0
  allocate(j3soil1(nx1,ny1,nzsoil1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyzsoil1,'j3soil1')
  j3soil1 = 1.0
  allocate(j3soilinv1(nx1,ny1,nzsoil1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyzsoil1,'j3soilinv1')
  j3soilinv1 = 1.0

  allocate(tem11(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'tem11')
  tem11 = 0.0
  allocate(tem21(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'tem21')
  tem21 = 0.0

  allocate(wgtsx(nx1,3),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*3,'wgtsx')
  wgtsx = 0.0
  allocate(wgtsy(ny1,3),stat=istatus)
  CALL alloc_status_accounting(istatus,ny1*3,'wgtsy')
  wgtsy = 0.0
  allocate(wgtux(nx1,3),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*3,'wgtux')
  wgtux = 0.0
  allocate(wgtvy(ny1,3),stat=istatus)
  CALL alloc_status_accounting(istatus,ny1*3,'wgtvy')
  wgtvy = 0.0
  allocate(wgtz(nx1,ny1,nz1,3),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1*3,'wgtz')
  wgtz  = 0.0

  allocate(isx(nx1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1,'isx')
  isx = 0.0
  allocate(jsy(ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,ny1,'jsy')
  jsy = 0.0
  allocate(iux(nx1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1,'iux')
  iux = 0.0
  allocate(jvy(ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,ny1,'jvy')
  jvy = 0.0
  allocate(kz(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxyz1,'kz')
  kz = 0.0
  allocate(zp1d1(nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nz1,'zp1d1')
  zp1d1 = 0.0
  allocate(dzp1d1(nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nz1,'dzp1d1')
  dzp1d1 = 0.0

  allocate(xs(nx ),stat=istatus)
  CALL alloc_status_accounting(istatus,nx,'xs')
  xs = 0.0
  allocate(ys(ny ),stat=istatus)
  CALL alloc_status_accounting(istatus,ny,'ys')
  ys = 0.0
  allocate(xs1(nx1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1,'xs1')
  xs1 = 0.0
  allocate(ys1(ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,ny1,'ys1')
  ys1 = 0.0

  allocate(xs_2d(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'xs_2d')
  xs_2d = 0.0
  allocate(ys_2d(nx,ny),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy,'ys_2d')
  ys_2d = 0.0

  allocate(temx1yz(nx1,ny ,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*ny*nz,'temx1yz')
  temx1yz = 0.0
  allocate(temx1y1z(nx1,ny1,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*ny1*nz,'temx1y1z')
  temx1y1z = 0.0
  allocate(temx1y1zb(nx1,ny1,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*ny1*nz,'temx1y1zb')
  temx1y1zb = 0.0

  allocate(temz1d1(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'temz1d1')
  temz1d1 = 0.0
  allocate(temz1d2(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'temz1d2')
  temz1d2 = 0.0
  allocate(temz1d3(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'temz1d3')
  temz1d3 = 0.0
  allocate(temz1d4(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'temz1d4')
  temz1d4 = 0.0
  allocate(temz1d5(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'temz1d5')
  temz1d5 = 0.0
  allocate(temz1d6(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'temz1d6')
  temz1d6 = 0.0
  allocate(temz1d7(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'temz1d7')
  temz1d7 = 0.0

  allocate(zsnd(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'zsnd')
  zsnd = 0.0
  allocate(ptsnd(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'ptsnd')
  ptsnd = 0.0
  allocate(qvsnd(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'qvsnd')
  qvsnd = 0.0
  allocate(psnd(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'psnd')
  psnd = 0.0
  allocate(ubsnd(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'ubsnd')
  ubsnd = 0.0
  allocate(vbsnd(lvlprof),stat=istatus)
  CALL alloc_status_accounting(istatus,lvlprof,'vbsnd')
  vbsnd = 0.0

  allocate(ptpsfc(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'ptpsfc')
  ptpsfc = 0.0
  allocate(ppsfc(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'ppsfc')
  ppsfc = 0.0
  allocate(qvsfc(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'qvsfc')
  qvsfc = 0.0
  allocate(qscalarsfc(nx1,ny1,nscalar),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1*nscalar,'qscalarsfc')
  qscalarsfc = 0.0
  allocate(tkesfc(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'tkesfc')
  tkesfc = 0.0
  allocate(kmhsfc(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'kmhsfc')
  kmhsfc = 0.0
  allocate(kmvsfc(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'kmvsfc')
  kmvsfc = 0.0
  allocate(ptbsfc(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'ptbsfc')
  ptbsfc = 0.0
  allocate(pbsfc(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'pbsfc')
  pbsfc = 0.0
  allocate(qvbsfc(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'qvbsfc')
  qvbsfc = 0.0
  allocate(htrnx1y1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'htrnx1y1')
  htrnx1y1 = 0.0
  allocate(radsfc(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'radsfc')
  radsfc = 0.0
  allocate(ktrnx1y1(nx1,ny1),stat=istatus)
  CALL alloc_status_accounting(istatus,nxy1,'ktrnx1y1')
  ktrnx1y1 = 0.0

  ALLOCATE (ix(nx1,ny1),stat=istatus)
  ALLOCATE (jy(nx1,ny1),stat=istatus)
  ALLOCATE (xw(nx1,ny1),stat=istatus)
  ALLOCATE (yw(nx1,ny1),stat=istatus)

  nxyz = nx*ny*nz

  u    (:,:,:) = 0.0
  v    (:,:,:) = 0.0
  uprt (:,:,:) = 0.0
  vprt (:,:,:) = 0.0
  w    (:,:,:) = 0.0
  ptprt(:,:,:) = 0.0
  pprt (:,:,:) = 0.0
  qv   (:,:,:) = 0.0
  qvprt(:,:,:) = 0.0
  CALL flzero(qscalar,nxyz*nscalar)
  tke  (:,:,:) = 0.0
  kmh  (:,:,:) = 0.0
  kmv  (:,:,:) = 0.0

  ubar  (:,:,:) = 0.0
  vbar  (:,:,:) = 0.0
  wbar  (:,:,:) = 0.0
  ptbar (:,:,:) = 0.0
  pbar  (:,:,:) = 0.0
  rhobar(:,:,:) = 0.0
  qvbar (:,:,:) = 0.0

  x(:) = 0.0
  y(:) = 0.0
  z(:) = 0.0
  zp(:,:,:) = 0.0

  tsoil  (:,:,:,:) = 0.0
  qsoil  (:,:,:,:) = 0.0
  wetcanp(:,:,:)   = 0.0

  soiltyp (:,:,:) = 0
  stypfrct(:,:,:) = 0.0

  vegtyp  = 0
  snowdpth= 0.0
  lai     = 0.0
  roufns  = 0.0
  veg     = 0.0
  ndvi    = 0.0
  raing   = 0.0
  rainc   = 0.0
  prcrate = 0.0
  prcrate = 0.0
  prcrate = 0.0
  prcrate = 0.0
  radsw   = 0.0
  radswnet = 0.0
  radlwin = 0.0
  rnflx   = 0.0
  usflx   = 0.0
  vsflx   = 0.0
  ptsflx  = 0.0
  qvsflx  = 0.0

  radfrc = 0.0

!
!-----------------------------------------------------------------------
!
!  Run interpolation.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Loop over data files
!
!-----------------------------------------------------------------------
!
  ireturn = 0

  DO nfile = 1,nhisfile

    filename = hisfile(nfile)

    lenfil=len_trim(filename)
    WRITE(6,'(/a,a,a)')                                                 &
        ' Data set ', filename(1:lenfil) ,' to be processed.'

    IF ( realtime /= 0 ) CALL check_file(filename,ntries,sleeptime)

!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
    102   CONTINUE

    CALL dtaread(nx,ny,nz,nzsoil,nstyp,                                 &
                 hinfmt, nchin,grdbasfn(1:lengbf),lengbf,               &
                 filename(1:lenfil),lenfil,time,                        &
                 x,y,z,zp,zpsoil,uprt,vprt,w,ptprt,pprt,                &
                 qvprt, qscalar, tke,kmh,kmv,                           &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 ireturn, tem1, tem2, tem3)

    IF (nstypin == 1) THEN
      stypfrct(:,:,1) = 1.
      IF (nstyp > nstypin) THEN
        DO is = 2,nstyp
          stypfrct(:,:,is) = 0.
        END DO
      END IF
    END IF

    curtim = time

    IF( hinfmt == 9 .AND. ireturn == 2 ) THEN
      WRITE(6,'(/1x,a/)') 'The end of GrADS file was reached.'
      CLOSE ( nchin )
      CALL retunit( nchin )
      GO TO 9001
    END IF

    IF( ireturn /= 0 ) GO TO 9002            ! Read was unsuccessful
!
!-----------------------------------------------------------------------
!
!  Grid and base related calculations that need to be done only once.
!
!-----------------------------------------------------------------------
!
    IF( nfile == 1) THEN ! Need to do this only once.
!
!-----------------------------------------------------------------------
!
!  Set-up land surface property data, if requested.
!
!-----------------------------------------------------------------------
!
      IF ( sfcdat == 0 ) GO TO 110

      IF ( sfcdat == 1 ) THEN

        nstyp = 1

        DO j=1, ny-1
          DO i=1, nx-1
            soiltyp(i,j,1) = styp
            vegtyp (i,j) = vtyp
            lai    (i,j) = lai0
            roufns (i,j) = roufns0
            veg    (i,j) = veg0
          END DO
        END DO

      ELSE IF (sfcdat == 2 .OR. (sfcdat == 3.AND.landin /= 1) ) THEN
!
!-----------------------------------------------------------------------
!
!  Read the surface property data from file sfcdtfl (passed
!  in globcst.inc).
!
!-----------------------------------------------------------------------
!

        ! Preserve history dump values.  We'll need to reconcile the
        ! different soil type contributions for tsoil, qsoil, wetcanp.

        soiltyp_his(:,:,:) = soiltyp(:,:,:)
        wetcanp_his(:,:,:) = wetcanp(:,:,:)
        tsoil_his(:,:,:,:) = tsoil(:,:,:,:)
        qsoil_his(:,:,:,:) = qsoil(:,:,:,:)

        stypfrct(:,:,:) = 0.0

        nstyp = nstyps_sfcdat

        CALL readsfcdt( nx,ny,nstyp,sfcdtfl,dx,dy,                      &
             mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,       &
             soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi )

        IF (nstyps_sfcdat == 1) THEN
          stypfrct(:,:,1) = 1.
          IF (nstyp > nstyps_sfcdat) THEN
            DO is = 2,nstyp
              stypfrct(:,:,is) = 0.
            END DO
          END IF
        END IF

        DO k = 1,nzsoil
          DO j = 1,ny-1
            DO i = 1,nx-1
              DO is = 1,nstyps_sfcdat
                tsoil(i,j,k,is) = tsoil(i,j,k,0)
                qsoil(i,j,k,is) = qsoil(i,j,k,0)
                IF (k == 1) wetcanp(i,j,is) = wetcanp(i,j,0)
                DO ii = 1,nstypin
                  IF (soiltyp_his(i,j,ii) == soiltyp(i,j,is)) THEN
                    tsoil(i,j,k,is) = tsoil_his(i,j,k,ii)
                    qsoil(i,j,k,is) = qsoil_his(i,j,k,ii)
                    IF (k == 1) wetcanp(i,j,is) = wetcanp_his(i,j,ii)
                  END IF
                END DO ! DO ii = 1,nstypin
              END DO ! DO is = 1,nstyps_sfcdat

              ! Now adjust the average values.
              tsoil(i,j,k,0) = 0.0
              qsoil(i,j,k,0) = 0.0
              IF (k == 1) wetcanp(i,j,0) = 0.0
              DO is = 1,nstyps_sfcdat
                tsoil(i,j,k,0) = tsoil(i,j,k,0) + &
                                 tsoil(i,j,k,is)*stypfrct(i,j,is)
                qsoil(i,j,k,0) = qsoil(i,j,k,0) + &
                                 qsoil(i,j,k,is)*stypfrct(i,j,is)
                IF (k == 1) wetcanp(i,j,0) = wetcanp(i,j,0) + &
                                 wetcanp(i,j,is)*stypfrct(i,j,is)
              END DO ! DO is = 1,nstyps_sfcdat

            END DO ! DO i = 1,nx-1
          END DO ! DO j = 1,ny-1
        END DO ! DO k = 1,nzsoil

        nstypin = nstyps_sfcdat

      ELSE IF(sfcdat == 3 .AND. landin == 1) THEN

        WRITE(6,'(1x,a/,a/)')                                           &
            'Surface property data in the history file was used.',      &
            'Data in ',sfcdtfl,' ignored.'

      ELSE

        WRITE(6,'(1x,a,i3,a/)')                                         &
            'Invalid surface data input option. sfcdat =',sfcdat,       &
            '. Program stopped in ARPSINTRP'
        STOP

      END IF  ! sfcdat option

      IF ( nstyp == 1 ) THEN
        DO j=1,ny
          DO i=1,nx
            stypfrct(i,j,1) = 1.0
          END DO
        END DO
      END IF

      110     CONTINUE

      DO j=1,ny
        DO i=1,nx
          hterain(i,j) = zp(i,j,2)
        END DO
      END DO

      dx = x(2) - x(1)
      dy = y(2) - y(1)

      xorig = x(2)
      yorig = y(2)
      zorig = z(2)

      WRITE(6,'(1x,a,3f15.3)') 'xorig, yorig, zorig =', xorig, yorig, zorig

      ebc = 3
      wbc = 3
      sbc = 3
      nbc = 3

      dx = x(2)-x(1)
      dy = y(2)-y(1)
      dz = z(2)-z(1)

    END IF
!
!-----------------------------------------------------------------------
!
!  Set up the map projection of the input grid.
!
!-----------------------------------------------------------------------
!
    alatpro(1)=trulat1
    alatpro(2)=trulat2

    dx = x(2)-x(1)
    dy = y(2)-y(1)
    IF( sclfct /= 1.0) THEN
      sclf  = 1.0/sclfct
      dxscl = dx*sclf
      dyscl = dy*sclf
    ELSE
      sclf  = 1.0
      dxscl = dx
      dyscl = dy
    END IF

    PRINT*,mapproj,sclf,alatpro,trulon,ctrlat,ctrlon

    IF (same_res == 1) THEN
      IF(dx /= dx1 .OR. dy /= dy1) THEN
        WRITE (*,*) "WARNING: resetting dx1 and/or dy1 to match input grid."
      ENDIF
      dx1 = dx
      dy1 = dy
    ENDIF

    ng = 0
500 CONTINUE ! Loop here for each output grid (noutgrds)
    IF (noutgrds > 0) THEN
      ng = ng + 1
      ctrlat1 = clat_grd(ng)
      ctrlon1 = clon_grd(ng)
      xctr1 = xctr_grd(ng)
      yctr1 = yctr_grd(ng)
    ENDIF

    CALL setmapr( mapproj,sclf,alatpro,trulon )
    CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
    swx = ctrx - (FLOAT(nx-3)/2.)*dxscl
    swy = ctry - (FLOAT(ny-3)/2.)*dyscl
    CALL setorig( 1, swx, swy) ! set up the model origin to the coord.

    CALL xytoll(1,1,(FLOAT(nx-3)/2.)*dx,(FLOAT(ny-3)/2.)*dy,            &
         tmp1,tmp2)
    PRINT*,'ctrlat=',tmp1,', ctrlon=',tmp2

    IF (xy_or_ll == 1) THEN
      CALL xytoll(1,1,xctr1,yctr1, ctrlat1,ctrlon1)
    ELSE IF (xy_or_ll == 2) THEN
      CALL lltoxy(1,1, ctrlat1,ctrlon1, xctr1,yctr1)
    ELSE
      WRITE (*,*) "ERROR: xy_or_ll =",xy_or_ll," not a valid value."
      STOP
    END IF


    IF (same_res == 1) THEN
      IF(dx /= dx1 .OR. dy /= dy1) THEN
        WRITE (*,*) "WARNING: resetting dx1 and/or dy1 to match input grid."
      ENDIF
      dx1 = dx
      dy1 = dy
    ENDIF

    ctrlat1_old = ctrlat1
    ctrlon1_old = ctrlon1

    IF( snap_to_grid /= 0 .OR. same_res == 1) THEN

      xsw = x(1)
      ysw = y(1)
      xsw1 = xctr1-(nx1-1)*dx1*0.5
      ysw1 = yctr1-(ny1-1)*dy1*0.5

      xsw1_old = xsw1
      ysw1_old = ysw1

      PRINT*,'mod =', MOD(dx,dx1),MOD(dy,dy1)

      IF(((dx >= dx1 .AND.                                              &
            ABS((MOD(dx,dx1)-dx1)*MOD(dx,dx1)) <= 0.001*dx1).OR.        &
            (dx <= dx1 .AND.                                            &
            ABS((MOD(dx1,dx)-dx )*MOD(dx1,dx)) <= 0.001*dx )).AND.      &
            ((dy >= dy1 .AND.                                           &
            ABS((MOD(dy,dy1)-dy1)*MOD(dy,dy1)) <= 0.001*dy1).OR.        &
            (dy <= dy1 .AND.                                            &
            ABS((MOD(dy1,dy)-dy )*MOD(dy1,dy)) <= 0.001*dy ) )) THEN
        PRINT*,'is multiple'

        IF(snap_to_grid == 1) THEN
          xsw1 = xsw + nint( (xsw1_old-xsw)/dx1 )*dx1
          ysw1 = ysw + nint( (ysw1_old-ysw)/dy1 )*dy1
        ELSE IF(snap_to_grid == 2) THEN
          xsw1 = xsw + nint( (xsw1_old-xsw)/dx )*dx
          ysw1 = ysw + nint( (ysw1_old-ysw)/dy )*dy
        END IF

        PRINT*,'xsw,ysw,xsw1,ysw1=',xsw,ysw,xsw1,ysw1
        PRINT*,'xsw1_old,ysw1_old=',xsw1_old,ysw1_old

        IF((snap_to_grid /= 0) .AND. (                                  &
              (ABS((xsw1_old - xsw1)) > 0.00001*dx1).OR.                &
              (ABS((ysw1_old - ysw1)) > 0.00001*dy1) )) THEN
          xctr1_old = xctr1
          yctr1_old = yctr1

          xctr1 = xsw1+(nx1-1)*dx1*0.5
          yctr1 = ysw1+(ny1-1)*dy1*0.5
!
!    Note there is a possibility that the shifted grid is out of bound
!    The checking is done later.
!

          WRITE(6,'(/1x,a/)')                                           &
              '################################################################'
          WRITE(6,'(/1x,a/,2(/1x,a)/)')                                 &
              '                    ATTENTION:',                         &
              'xctr1 and/or yctr1 reset so that some of the output ',   &
              'grid lines match those of input grid.'

          WRITE(6,'(1x,2(a,f17.5))')'Original xctr1_old=',xctr1_old,    &
                                  ', yctr1_old=',yctr1_old
          WRITE(6,'(1x,2(a,f17.5))')'New      xctr1   = ',xctr1,        &
                                  ', yctr1   = ',yctr1
          WRITE(6,'(/1x,a/)')                                           &
              '################################################################'

          CALL xytoll(1,1,xctr1,yctr1, ctrlat1,ctrlon1)
        END IF
      ELSE
        same_res = 0  ! ???
      END IF

    END IF

    IF( nz1 /= nz .AND. same_res == 1) THEN
      WRITE (6,'(/,1x,a,I3,a,I3,a,/,10x,a,/)')                          &
          'WARNING: nz1 (',nz1,') does not match nz (',nz,').',         &
          'Resetting same_res to 0.'
      same_res = 0
    ENDIF
    IF( nzsoil1 /= nzsoil .AND. same_res == 1) THEN
      WRITE (*,*) "WARNING: nzsoil1 (",nzsoil1,") does not match nzsoil (",nzsoil,")."
      WRITE (*,*) "   Resetting same_res to 0."
      same_res = 0
    ENDIF

    IF (same_res == 1) THEN
      IF(dzsoil /= dzsoil1 ) THEN
        WRITE (*,*) "WARNING: resetting dzsoil1 to match input grid."
      ENDIF
      dzsoil1 = dzsoil
    ENDIF


    WRITE(6,'(1x,2(a,f17.5))')                                          &
        'Original ctrlat1_old = ',ctrlat1_old,', ctrlon1_old = ',       &
        ctrlon1_old
    WRITE(6,'(1x,2(a,f17.5))')                                          &
        'Final ctrlat1        = ',ctrlat1    ,', ctrlon1     = ',ctrlon1
    WRITE(6,'(1x,2(a,f17.5))')                                          &
        'Final xctr1          = ',xctr1,      ', yctr1       = ',yctr1
    WRITE(6,'(1x,a,i2/)') 'same_res = ', same_res

!
!-----------------------------------------------------------------------
!
!  Perform spatial interpolations.
!
!-----------------------------------------------------------------------
!
!    IF( nfile == 1) THEN
    IF( nfile == 1 .or. noutgrds > 0) THEN
       ! Need to do this for the first file only

      xorig1 = xctr1 - (nx1-3)*dx1*0.5
      yorig1 = yctr1 - (ny1-3)*dy1*0.5
      zorig1 = zorig

      DO i=1,nx1
        x1(i) =xorig1+(i-2.0)*dx1
        xs1(i)=xorig1+(i-1.5)*dx1
      END DO

      DO j=1,ny1
        y1 (j) =yorig1+(j-2.0)*dy1
        ys1(j) =yorig1+(j-1.5)*dy1
      END DO

      ib= nint(x1(1)/dx)+1
      jb= nint(y1(1)/dy)+1

      DO i=1,nx-1
        xs(i)=0.5*(x(i)+x(i+1))
      END DO

      DO j=1,ny-1
        ys(j)=0.5*(y(j)+y(j+1))
      END DO

      DO k=1,nz1
        z1(k)=zorig1+(k-2.0)*dz1
      END DO

      IF( nfile == 1) THEN  ! Need to do this for the first file only
        xeps = 0.01*dx
        yeps = 0.01*dy

        WRITE(6,'(1x,a,4f15.2)') 'x (2),x (nx -1),y (2),y (ny -1) =',     &
                x(2),x(nx-1),y(2),y(ny-1)
        WRITE(6,'(1x,a,4F15.2)') 'x1(2),x1(nx1-1),y1(2),y1(ny1-1) =',     &
                x1(2),x1(nx1-1),y1(2),y1(ny1-1)

        IF(x1(    1) < x(   1)-xeps.OR.x1(nx1) > x(nx)+xeps.OR.           &
              y1(    1) < y(   1)-yeps.OR.y1(ny1) > y(ny)+yeps) THEN
          WRITE(6,'(3(/2x,a),/2x,2(a,f12.4),2(a,i5),2(a,f12.4),/2x,a)')   &
        'Sorry, at least part of your new grid is outside the border of', &
              'the original grid, please check input parameters',         &
              'dx1,dy1, nx1, ny1, xctr1 and yctr1. Currently,',           &
              'dx1=',dx1,', dy1=',dy1,', nx1=',nx1,', ny1=',ny1,          &
              ', xctr1=',xctr1,', yctr1=',yctr1,                          &
              'Job stopped in ARPSINTRP.'
          GOTO 600
!          STOP 1001
        END IF
      END IF

      IF( same_res == 0 ) THEN

        IF( intrphopt == 1) THEN ! Weights for 1st order interpolation

          DO i=1,nx1-1
            isx(i) = MAX(1, MIN(nx-2, INT((xs1(i)-xs(1))/dx)+1 ))
            wgtsx(i,1)= (xs(isx(i)+1)-xs1(i))/(xs(isx(i)+1)-xs(isx(i)))
          END DO

          DO j=1,ny1-1
            jsy(j) = MAX(1, MIN(ny-2, INT((ys1(j)-ys(1))/dy)+1 ))
            wgtsy(j,1)= (ys(jsy(j)+1)-ys1(j))/(ys(jsy(j)+1)-ys(jsy(j)))
          END DO

          DO i=1,nx1
            iux(i) = MAX(1, MIN(nx-1, INT((x1 (i)-x (1))/dx)+1 ))
            wgtux(i,1)= (x (iux(i)+1)-x1 (i))/(x (iux(i)+1)-x (iux(i)))
          END DO

          DO j=1,ny1
            jvy(j) = MAX(1, MIN(ny-1, INT((y1 (j)-y (1))/dy)+1 ))
            wgtvy(j,1)= (y (jvy(j)+1)-y1 (j))/(y (jvy(j)+1)-y (jvy(j)))
          END DO

        ELSE ! Weights for second order interpolation

          DO i=1,nx1-1
            isx(i) = MAX(2, MIN(nx-2, NINT((xs1(i)-xs(1))/dx)+1 ))
            d=xs1(i)
            a=xs(isx(i)-1)
            b=xs(isx(i))
            c=xs(isx(i)+1)
            wgtsx(i,1)= (d-b)*(d-c)/((a-b)*(a-c))
            wgtsx(i,2)= (d-a)*(d-c)/((b-a)*(b-c))
            wgtsx(i,3)= (d-a)*(d-b)/((c-a)*(c-b))
          END DO

          DO j=1,ny1-1
            jsy(j) = MAX(2, MIN(ny-2, NINT((ys1(j)-ys(1))/dy)+1 ))
            d=ys1(j)
            a=ys(jsy(j)-1)
            b=ys(jsy(j))
            c=ys(jsy(j)+1)
            wgtsy(j,1)= (d-b)*(d-c)/((a-b)*(a-c))
            wgtsy(j,2)= (d-a)*(d-c)/((b-a)*(b-c))
            wgtsy(j,3)= (d-a)*(d-b)/((c-a)*(c-b))
          END DO

          DO i=1,nx1
            iux(i) = MAX(2, MIN(nx-1, NINT((x1 (i)-x (1))/dx)+1 ))
            d=x1(i)
            a=x(iux(i)-1)
            b=x(iux(i))
            c=x(iux(i)+1)
            wgtux(i,1)= (d-b)*(d-c)/((a-b)*(a-c))
            wgtux(i,2)= (d-a)*(d-c)/((b-a)*(b-c))
            wgtux(i,3)= (d-a)*(d-b)/((c-a)*(c-b))
          END DO

          DO j=1,ny1
            jvy(j) = MAX(2, MIN(ny-1, NINT((y1 (j)-y (1))/dy)+1 ))
            d=y1(j)
            a=y(jvy(j)-1)
            b=y(jvy(j))
            c=y(jvy(j)+1)
            wgtvy(j,1)= (d-b)*(d-c)/((a-b)*(a-c))
            wgtvy(j,2)= (d-a)*(d-c)/((b-a)*(b-c))
            wgtvy(j,3)= (d-a)*(d-b)/((c-a)*(c-b))
          END DO

        END IF ! intrphopt = ?


        IF( ternopt1 == 0 ) THEN

          DO i=1,nx1-1
            DO j=1,ny1-1
              hterain1(i,j) = zrefsfc1
            END DO
          END DO

        ELSE IF( ternopt1 == 1 ) THEN  ! Bell-shaped mountain
!
!-----------------------------------------------------------------------
!
!  Define the bell-shaped mountain
!
!-----------------------------------------------------------------------
!
          pi2 = 1.5707963267949

          DO j=1,ny1-1
            DO i=1,nx1-1

              IF( mntwidy1 < 0.0 ) THEN   ! 2-d terrain in x-z plane.

                radnd = 1.0+((xs1(i)-mntctrx1)/mntwidx1)**2

              ELSE IF( mntwidx1 < 0.0 ) THEN ! 2-d terrain in y-z plane.

                radnd = 1.0+((ys1(j)-mntctry1)/mntwidy1)**2

              ELSE                             ! 3-d terrain

                radnd = 1.0+((xs1(i)-mntctrx1)/mntwidx1)**2             &
                       +((ys1(j)-mntctry1)/mntwidy1)**2
              END IF

              hterain1(i,j) = zrefsfc1 + hmount1/radnd

            END DO
          END DO

        ELSE IF( ternopt1 == 2 .or. ternopt1 == 4 ) THEN
                                              ! Read from terrain data base

!
!-----------------------------------------------------------------------
!
!    Read in the terrain data.
!
!-----------------------------------------------------------------------
!
          lenstr = 256
          CALL strlnth( terndta1, lenstr)
          CALL readtrn( nx1,ny1,dx1,dy1, terndta1(1:lenstr),            &
                 mapproj,trulat1,trulat2,trulon,sclfct,ctrlat1,ctrlon1, &
                 hterain1 )

!
!-----------------------------------------------------------------------
!
!  Define the new grid terrain by interpolation
!
!-----------------------------------------------------------------------
!
        END IF

        IF( ternopt1 == 3 .or. ternopt1 == 4 ) THEN

          CALL intrpxy3d(hterain,nx,1,nx-1,ny,1,ny-1,1,1,1,             &
                         wgtsx,isx,wgtsy,jsy,intrphopt,                 &
                         htrn1orig,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

          !hterain1 = htrn1orig   ! uncomment if ternopt1=4 (ntmerge) dissabled
          IF (ternopt1 == 3) THEN
            hterain1 = htrn1orig
          ELSE
            ntmergeinv = 1.d0/ntmerge
            DO j = 1,ny1
              DO i = 1,nx1
                idist = max(0,min(ntmerge,i-2,nx1-2-i,j-2,ny1-2-j))
                mfac = idist*ntmergeinv
                hterain1(i,j) = (1.d0-mfac)*htrn1orig(i,j)  &
                                + mfac*hterain1(i,j)
              END DO
            END DO
          END IF

        END IF

        npoint_below_ground = 0

        IF( ternopt1 /= 3) THEN
          CALL intrpxy3d(hterain,nx,1,nx-1,ny,1,ny-1,1,1,1,             &
                         wgtsx,isx,wgtsy,jsy,intrphopt,                 &
                         temx1y1z,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

          CALL a3dmax0(hterain,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,         &
                      amax,amin)
          WRITE(6,'(1x,2(a,e13.6))') 'htmin   = ', amin,', htmax    =',amax

          CALL a3dmax0(temx1y1z,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,    &
                      amax,amin)
          WRITE(6,'(1x,2(a,e13.6))') 'ht1min  = ', amin,', ht1max   =',amax

          CALL a3dmax0(hterain1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,    &
                      amax,amin)
          WRITE(6,'(1x,2(a,e13.6))') 'ht1min  = ', amin,', ht1max   =',amax

          npoint_below_ground = 0
          DO i=1,nx1-1
            DO j=1,ny1-1
              IF( hterain1(i,j) < temx1y1z(i,j,1)-1.0E-2) THEN
                npoint_below_ground = 1
                GO TO 2100
              END IF
            END DO
          END DO
          2100    CONTINUE
        END IF

!
!-----------------------------------------------------------------------
!
!  Set up a stretched vertical grid for the output grid.
!
!  For strhopt1=1, function y = a+b*x**3 is used to specify dz as a
!                              function of k.
!  For strhopt1=2, function y = c + a*tanh(b*x) is used to specify dz
!                              as a function of k.
!
!-----------------------------------------------------------------------
!

        IF ( strhopt1 == 0 ) THEN

          DO k=1,nz1
            zp1d1(k) = z1(k)
          END DO

        ELSE IF ( strhopt1 == 1 .OR.strhopt1 == 2 ) THEN

          za = zrefsfc1 + MAX(0.0, MIN(dlayer11, z1(nz1-2)-zrefsfc1 ))
          zb = za       + MAX(0.0, MIN(dlayer21, z1(nz1-1)-za      ))

          IF( dlayer11 >= (nz1-3)*dzmin1 ) THEN
            WRITE(6,'(/1x,a,f13.3,/a,f13.3,a,a)')                       &
                'Can not setup a vertical grid with uniform dz=',dzmin1, &
                ' over the depth of ',dlayer11,' please specify a smaller ', &
                'value of dlayer11. Program stopped in ARPSINTRP.'
            STOP 9105
          END IF

          CALL strhgrd(nz1,strhopt1,zrefsfc1,za,zb,z1(nz1-1),           &
                       dzmin1,strhtune1, zp1d1,dzp1d1)

        ELSE

          WRITE(6,'(1x,a,i3,a/)')                                       &
              'Invalid vertical grid stretching option, strhopt1 was ',strhopt1, &
              '. Program stopped in INIGRD.'
          STOP 9106

        END IF

!-----------------------------------------------------------------------
!
!  Physical height of computational grid defined as
!
!  Zp=(z-zrefsfc)*(Zm-hterain)/(Zm-zrefsfc)+hterain for z=<Zm.
!  ZP=z for z>Zm
!
!  where Zm the height at which the grid levels becomes flat.
!  Hm < Zm =< Ztop, hm is the height of mountain and Ztop the height
!  of model top.
!
!-----------------------------------------------------------------------
!
        DO k=nz1-1,2,-1
          IF(zp1d1(k) <= zflat1) THEN
            zflat11 = zp1d1(k)
            EXIT
          END IF
        END DO
!        525   CONTINUE
        zflat11=MAX(MIN(z1(nz1-1),zflat11),zrefsfc1)

        DO k=2,nz1-1

          IF(zp1d1(k) > zflat11) THEN
            DO j=1,ny1-1
              DO i=1,nx1-1
                zp1(i,j,k)=zp1d1(k)
              END DO
            END DO
          ELSE
            DO j=1,ny1-1
              DO i=1,nx1-1
                zp1(i,j,k)=(zp1d1(k)-zrefsfc1)*(zflat11-hterain1(i,j))  &
                           /(zflat11-zrefsfc1)+hterain1(i,j)
              END DO
            END DO
          END IF

        END DO

        DO j=1,ny1-1
          DO i=1,nx1-1
            zp1(i,j,2)=hterain1(i,j)
            zp1(i,j,1)=2.0*zp1(i,j,2)-zp1(i,j,3)
            zp1(i,j,nz1)=2.0*zp1(i,j,nz1-1)-zp1(i,j,nz1-2)
          END DO
        END DO
!
        CALL jacob(nx1,ny1,nz1,x1,y1,z1,zp1,j11,j21,j31,tem1)

        dzsoil_tmp = dzsoil
        dzsoil = dzsoil1 ! Copy to global variable
        CALL inisoilgrd(nx1,ny1,nzsoil1,hterain1,zpsoil1, &
                        j3soil1,j3soilinv1)
        dzsoil = dzsoil_tmp ! Restore original value

      END IF ! same_res = 0?

    END IF  ! IF (nfile.eq.1)

!
!-----------------------------------------------------------------------
!
!  Calculate total fields from that for base state and perturbations
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          u (i,j,k)=uprt (i,j,k)+ubar (i,j,k)
          v (i,j,k)=vprt (i,j,k)+vbar (i,j,k)
          qv(i,j,k)=qvprt(i,j,k)+qvbar(i,j,k)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
!  If the grids are of the same resolution, begin copying variables
!  from external grid to new grid.
!
!-----------------------------------------------------------------------

    IF( same_res == 1 ) THEN

      DO k=1,nz1
        DO j=1,ny1
          DO i=1,nx1
            u1    (i,j,k)=u    (i+ib,j+jb,k)
            v1    (i,j,k)=v    (i+ib,j+jb,k)
            w1    (i,j,k)=w    (i+ib,j+jb,k)
            ptprt1(i,j,k)=ptprt(i+ib,j+jb,k)
            pprt1 (i,j,k)=pprt (i+ib,j+jb,k)
            qv1   (i,j,k)=qv   (i+ib,j+jb,k)
            DO nq=1,nscalar
              qscalar1(i,j,k,nq)=qscalar(i+ib,j+jb,k,nq)
            END DO
            tke1  (i,j,k)=tke  (i+ib,j+jb,k)
            kmh1  (i,j,k)=kmh  (i+ib,j+jb,k)
            kmv1  (i,j,k)=kmv  (i+ib,j+jb,k)
            radfrc1(i,j,k)=radfrc(i+ib,j+jb,k)
            rhobar1(i,j,k)=rhobar(i+ib,j+jb,k)
          END DO
        END DO
      END DO

!      IF( nfile == 1) THEN
      IF( nfile == 1 .or. noutgrds > 0) THEN

        DO k=1,nz1
          DO i=1,nx1
            DO j=1,ny1
              ubar1 (i,j,k)=ubar (i+ib,j+jb,k)
              vbar1 (i,j,k)=vbar (i+ib,j+jb,k)
              ptbar1(i,j,k)=ptbar(i+ib,j+jb,k)
              pbar1 (i,j,k)=pbar (i+ib,j+jb,k)
              qvbar1(i,j,k)=qvbar(i+ib,j+jb,k)
              zp1   (i,j,k)=zp   (i+ib,j+jb,k)
            END DO
          END DO
        END DO

        CALL jacob(nx1,ny1,nz1,x1,y1,z1,zp1,j11,j21,j31,tem1)

        dzsoil_tmp = dzsoil
        dzsoil = dzsoil1         ! Copy to global variable
        CALL inisoilgrd(nx1,ny1,nzsoil1,hterain1,zpsoil1, &
                        j3soil1,j3soilinv1)
        dzsoil = dzsoil_tmp      ! Restore original value

      END IF

      WRITE(*,*)'NOTE:  Ignoring newgrid_soil settings (same_res = 1).'

      DO i=1,nx1
        DO j=1,ny1
          hterain1(i,j)=hterain(i+ib,j+jb)

          IF (nstyp1 == nstypin) THEN ! Same soil types
            soiltyp1 (i,j,:)=soiltyp (i+ib,j+jb,:)
            stypfrct1(i,j,:)=stypfrct(i+ib,j+jb,:)
            wetcanp1 (i,j,:)=wetcanp (i+ib,j+jb,:)
            tsoil1   (i,j,:,:)=tsoil   (i+ib,j+jb,:,:)
            qsoil1   (i,j,:,:)=qsoil   (i+ib,j+jb,:,:)
          ELSE IF (nstyp1 > nstypin) THEN ! Expect more soil types
                                          ! than available.
            DO is = 1,nstypin
              soiltyp1 (i,j,is)=soiltyp (i+ib,j+jb,is)
              stypfrct1(i,j,is)=stypfrct(i+ib,j+jb,is)
              wetcanp1 (i,j,is)=wetcanp (i+ib,j+jb,is)
              tsoil1   (i,j,:,is)=tsoil   (i+ib,j+jb,:,is)
              qsoil1   (i,j,:,is)=qsoil   (i+ib,j+jb,:,is)
            END DO
            tsoil1(i,j,:,0) = tsoil   (i+ib,j+jb,:,0)
            qsoil1(i,j,:,0) = qsoil   (i+ib,j+jb,:,0)
            wetcanp1(i,j,0) = wetcanp (i+ib,j+jb,0)
            DO is = nstypin+1,nstyp1
              soiltyp1 (i,j,is)=soiltyp1(i,j,is-1)
              stypfrct1(i,j,is)=0.
              wetcanp1 (i,j,is)=wetcanp1(i,j,0)
              tsoil1   (i,j,:,is)=tsoil1(i,j,:,0)
              qsoil1   (i,j,:,is)=qsoil1(i,j,:,0)
            END DO
            tsoil1(i,j,:,0) = 0.
            qsoil1(i,j,:,0) = 0.
            wetcanp1(i,j,0) = 0.
            DO is = 1,nstyp1
              wetcanp1(i,j,0) = wetcanp1(i,j,0) + &
                                stypfrct1(i,j,is)*wetcanp1(i,j,is)
              tsoil1(i,j,:,0) = tsoil1(i,j,:,0) + &
                                stypfrct1(i,j,is)*tsoil1(i,j,:,is)
              qsoil1(i,j,:,0) = qsoil1(i,j,:,0) + &
                                stypfrct1(i,j,is)*qsoil1(i,j,:,is)
            END DO ! DO is = 1,nstyp1
          ELSE ! nstyp1 < nstypin         ! Use subset of available
                                          ! soil types.
            frctot = 0.
            DO is = 1,nstyp1
              soiltyp1 (i,j,is)=soiltyp (i+ib,j+jb,is)
              stypfrct1(i,j,is)=stypfrct(i+ib,j+jb,is)
              wetcanp1 (i,j,is)=wetcanp (i+ib,j+jb,is)
              tsoil1   (i,j,:,is)=tsoil   (i+ib,j+jb,:,is)
              qsoil1   (i,j,:,is)=qsoil   (i+ib,j+jb,:,is)
              frctot = frctot + stypfrct1(i,j,is)
            END DO
            IF (frctot /= 0) THEN
              DO is = 1,nstyp1
                stypfrct1(i,j,is) = stypfrct1(i,j,is)/frctot
              END DO ! DO is = 1,nstyp1
            ELSE
              stypfrct1(i,j,1) = 1.
              IF (nstyp1 > 1) THEN
                DO is = 2,nstyp1
                  stypfrct1(i,j,is) = 0.
                END DO
              END IF
            END IF
            tsoil1(i,j,:,0) = 0.
            qsoil1(i,j,:,0) = 0.
            wetcanp1(i,j,0) = 0.
            DO is = 1,nstyp1
              wetcanp1(i,j,0) = wetcanp1(i,j,0) + &
                                stypfrct1(i,j,is)*wetcanp1(i,j,is)
              tsoil1(i,j,:,0) = tsoil1(i,j,:,0) + &
                                stypfrct1(i,j,is)*tsoil1(i,j,:,is)
              qsoil1(i,j,:,0) = qsoil1(i,j,:,0) + &
                                stypfrct1(i,j,is)*qsoil1(i,j,:,is)
            END DO ! DO is = 1,nstyp1
          END IF

          vegtyp1(i,j)=vegtyp(i+ib,j+jb)
          lai1   (i,j)=lai   (i+ib,j+jb)
          roufns1(i,j)=roufns(i+ib,j+jb)
          veg1   (i,j)=veg   (i+ib,j+jb)

          snowdpth1(i,j)=snowdpth(i+ib,j+jb)

          raing1  (i,j)=raing(i+ib,j+jb)
          rainc1  (i,j)=rainc(i+ib,j+jb)

          prcrate1(i,j,1:4)=prcrate(i+ib,j+jb,1:4)

          radsw1 (i,j)=radsw (i+ib,j+jb)
          rnflx1 (i,j)=rnflx (i+ib,j+jb)
          radswnet1 (i,j)=radswnet(i+ib,j+jb)
          radlwin1 (i,j)=radlwin(i+ib,j+jb)
          usflx1 (i,j)=usflx (i+ib,j+jb)
          vsflx1 (i,j)=vsflx (i+ib,j+jb)
          ptsflx1(i,j)=ptsflx(i+ib,j+jb)
          qvsflx1(i,j)=qvsflx(i+ib,j+jb)
        END DO
      END DO

      GO TO 800 ! Skip interpolations

    END IF
!
!-----------------------------------------------------------------------
!
!   If different resolutions, interpolate the 2-D arrays to the new
!   grid.
!
!-----------------------------------------------------------------------
!

    CALL intrpxy3d(rainc,nx,1,nx-1,ny,1,ny-1,1,1,1,                     &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   rainc1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(raing,nx,1,nx-1,ny,1,ny-1,1,1,1,                     &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   raing1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(prcrate,nx,1,nx-1,ny,1,ny-1,4,1,4,                   &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   prcrate1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(radsw,nx,1,nx-1,ny,1,ny-1,1,1,1,                     &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   radsw1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(rnflx,nx,1,nx-1,ny,1,ny-1,1,1,1,                     &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   rnflx1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(radswnet,nx,1,nx-1,ny,1,ny-1,1,1,1,                  &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   radswnet1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(radlwin,nx,1,nx-1,ny,1,ny-1,1,1,1,                   &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   radlwin1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(usflx,nx,1,nx-1,ny,1,ny-1,1,1,1,                     &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   usflx1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(vsflx,nx,1,nx-1,ny,1,ny-1,1,1,1,                     &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   vsflx1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(ptsflx,nx,1,nx-1,ny,1,ny-1,1,1,1,                    &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   ptsflx1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(qvsflx,nx,1,nx-1,ny,1,ny-1,1,1,1,                    &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   qvsflx1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(lai,nx,1,nx-1,ny,1,ny-1,1,1,1,                       &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   lai1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(roufns,nx,1,nx-1,ny,1,ny-1,1,1,1,                    &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   roufns1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    CALL intrpxy3d(veg,nx,1,nx-1,ny,1,ny-1,1,1,1,                       &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   veg1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

!   Make sure external data and input nzsoil1 will
!   permit the use of the old ARPS Force-Restore Soil Model.

    IF (soilmodel_option == 1) THEN
      IF (nzsoil /= 2) THEN
        WRITE(6,*)"WARNING:  nzsoil from external ARPS grid is not 2"
        WRITE(6,*)"Will reset soilmodel_option to 2 (new soil model)"
        soilmodel_option = 2
      ELSE IF (nzsoil1 /= 2) THEN
        WRITE(6,*)"WARNING:  nzsoil1 from input file is not 2"
        WRITE(6,*)"Will reset soilmodel_option to 2 (new soil model)"
        soilmodel_option = 2
      END IF
    END IF

    IF (soilmodel_option == 1) THEN ! Old ARPS Force-Restore Soil Model

      ALLOCATE(tsoiltem1(nx,ny,0:nstyp),stat=istatus)
      CALL alloc_status_accounting(istatus,nxy*(nstyp+1),'tsoiltem1  ')
      ALLOCATE(tsoiltem2(nx,ny,0:nstyp),stat=istatus)
      CALL alloc_status_accounting(istatus,nxy*(nstyp+1),'tsoiltem2  ')
      ALLOCATE(qsoiltem1(nx,ny,0:nstyp),stat=istatus)
      CALL alloc_status_accounting(istatus,nxy*(nstyp+1),'qsoiltem1  ')
      ALLOCATE(qsoiltem2(nx,ny,0:nstyp),stat=istatus)
      CALL alloc_status_accounting(istatus,nxy*(nstyp+1),'qsoiltem2  ')

      ALLOCATE(tsoil1tem1(nx1,ny1,0:nstyp1),stat=istatus)
      CALL alloc_status_accounting(istatus,nxy1*(nstyp1+1),'tsoil1tem1  ')
      ALLOCATE(tsoil1tem2(nx1,ny1,0:nstyp1),stat=istatus)
      CALL alloc_status_accounting(istatus,nxy1*(nstyp1+1),'tsoil1tem2  ')
      ALLOCATE(qsoil1tem1(nx1,ny1,0:nstyp1),stat=istatus)
      CALL alloc_status_accounting(istatus,nxy1*(nstyp1+1),'qsoil1tem1  ')
      ALLOCATE(qsoil1tem2(nx1,ny1,0:nstyp1),stat=istatus)
      CALL alloc_status_accounting(istatus,nxy1*(nstyp1+1),'qsoil1tem2  ')

      tsoiltem1(:,:,:) = tsoil(:,:,1,:)
      tsoiltem2(:,:,:) = tsoil(:,:,2,:)
      qsoiltem1(:,:,:) = qsoil(:,:,1,:)
      qsoiltem2(:,:,:) = qsoil(:,:,2,:)

      DO j = 1,ny1-1
        DO i = 1,nx1-1
          ix(i,j) = isx(i)
          jy(i,j) = jsy(j)

          xw(i,j) =  &
           MAX(0.,MIN(1.,(xs(isx(i)+1)-xs1(i))/(xs(isx(i)+1)-xs(isx(i)))))
          yw(i,j) =  &
           MAX(0.,MIN(1.,(ys(jsy(j)+1)-ys1(j))/(ys(jsy(j)+1)-ys(jsy(j)))))
        END DO
      END DO

      WRITE(6,*) 'Calling intrp_soil ...'
      CALL intrp_soil(nx,ny,nx1,ny1,nstypin,nstyp1,xw,yw,ix,jy,           &
                    tsoiltem1,tsoiltem2,qsoiltem1,qsoiltem2, wetcanp,     &
                    soiltyp,stypfrct,vegtyp,                              &
                    tsoil1tem1,tsoil1tem2,qsoil1tem1,qsoil1tem2,wetcanp1, &
                    soiltyp1,stypfrct1,vegtyp1)

      tsoil1(:,:,1,:) = tsoil1tem1(:,:,:)
      tsoil1(:,:,2,:) = tsoil1tem2(:,:,:)
      qsoil1(:,:,1,:) = qsoil1tem1(:,:,:)
      qsoil1(:,:,2,:) = qsoil1tem2(:,:,:)

      DEALLOCATE (tsoiltem1, tsoiltem2)
      DEALLOCATE (qsoiltem1, qsoiltem2)
      DEALLOCATE (tsoil1tem1, tsoil1tem2)
      DEALLOCATE (qsoil1tem1, qsoil1tem2)

    ELSE                                              ! New OUSoil Model

      ! Convert to soil depth
      DO k = 1,nzsoil
        DO j =1,ny
          DO i =1,nx
            zpsoil(i,j,k) = hterain(i,j) - zpsoil(i,j,k)
          END DO
        END DO
      END DO
      DO k = 1,nzsoil1
        DO j =1,ny1
          DO i =1,nx1
            zpsoil1(i,j,k) = hterain1(i,j) - zpsoil1(i,j,k)
          END DO
        END DO
      END DO

      DO j=1,ny1
        DO i=1,nx1
          x2d1(i,j) = xs1(i)
          y2d1(i,j) = ys1(j)
        END DO ! DO i=1,nx1
      END DO ! DO j=1,ny1

      CALL setdxdy(nx,ny,1,nx,1,ny,x,y,dxfld,dyfld,rdxfld,rdyfld)

      DO k = 1,nzsoil-1
        DO j =1,ny-1
          DO i =1,nx-1
            rdzsoilfld(i,j,k) = &
              zpsoil(i,j,k+1) - zpsoil(i,j,k)
          END DO
        END DO
      END DO

      CALL intrpsoil3d_pst(nx,ny,nzsoil,nstypin,soiltyp,stypfrct, &
                           vegtyp,tsoil,qsoil,wetcanp,x,y,zpsoil, &
                           rdxfld,rdyfld,rdzsoilfld, &
                           nx1,ny1,nzsoil1,nstyp1,soiltyp1,stypfrct1, &
                           vegtyp1,tsoil1,qsoil1,wetcanp1,x2d1,y2d1, &
                           zpsoil1,i2d1,j2d1,k3d1)

      ! Convert back to MSL m
      DO k = 1,nzsoil
        DO j =1,ny
          DO i =1,nx
            zpsoil(i,j,k) = hterain(i,j) - zpsoil(i,j,k)
          END DO
        END DO
      END DO
      DO k = 1,nzsoil1
        DO j =1,ny1
          DO i =1,nx1
            zpsoil1(i,j,k) = hterain1(i,j) - zpsoil1(i,j,k)
          END DO
        END DO
      END DO

    END IF ! Force-restore or OU Soil.

    CALL intrpxy3d(snowdpth,nx,1,nx-1,ny,1,ny-1,1,1,1,                  &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   snowdpth1,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

!
!    DO i=1,nx1-1
!      DO j=1,ny1-1
!        vegtyp1 (i,j) = vegtyp (isx(i),jsy(j))
!      END DO
!    END DO
!
!    DO iss=1,nstyp
!      DO j=1,ny1-1
!        DO i=1,nx1-1
!          soiltyp1 (i,j,iss) = soiltyp (isx(i),jsy(j),iss)
!          stypfrct1(i,j,iss) = stypfrct(isx(i),jsy(j),iss)
!        END DO
!      END DO
!    END DO

!
!-----------------------------------------------------------------------
!
!   Interpolate 3D arrays
!
!-----------------------------------------------------------------------
!
    IF( z_or_p == 2) THEN  ! This situation has not been tested since
                           ! last major change to this program

!      ALLOCATE(zp_saved (nx, ny, nz),  STAT = istatus)
!      ALLOCATE(zp1_saved(nx1,ny1,nz1), STAT = istatus)
!
!      zp_saved(:,:,:) = zp(:,:,:)

      DO k=1,nz
        DO j=1,ny-1
          DO i=1,nx-1
            zp(i,j,k) = -ALOG( (pbar(i,j,k)+pprt(i,j,k))*1.0E-5 )
          END DO
        END DO
      END DO

      ptop = -ALOG(plevels(nz1)*1.0E-3)
      IF ( ANY(zp(:,:,nz-1) < ptop) ) THEN
        WRITE(6,'(/,1x,a,f7.1,a,/,8x,a,f7.1,a)')                        &
        'ERROR: The highest pressure level (',plevels(nz1),             &
        ') exceed the original model','top pressure at some location (',&
        EXP(-1*MINVAL(zp(:,:,nz-1)))*1.0E3,').'
        WRITE(6,'(8x,a,/)') 'Please change parameter "plevels" and try again.'
        STOP
      END IF

      DO k=1,nz1
        DO j=1,ny1-1
          DO i=1,nx1-1
            zp1(i,j,k) = -ALOG( plevels(k)*1.0E-3 )
          END DO
        END DO
      END DO
    END IF


    IF( npoint_below_ground /= 0 .AND. z_or_p == 1 .AND.                &
        (bglopt == 4 .OR. bglopt == 5) ) THEN
          ! Reconstruct horizontally homogeneous base state
      redo_base_state = 1
    ELSE
      redo_base_state = 0
    END IF

    IF( redo_base_state /= 0 ) THEN

      zpmin=0.5*(zp(2,2, 2) + zp(2,2,3)   )
      zpmax=0.5*(zp(2,2,nz) + zp(2,2,nz-1))
      DO i=2,nx-2
        DO j=2,ny-2
          zpmin=MIN(zpmin,0.5*(zp(i,j, 2)+zp(i,j,   3)))
          zpmax=MAX(zpmax,0.5*(zp(i,j,nz)+zp(i,j,nz-1)))
        END DO
      END DO

      zpmin1=MIN(zpmin1,0.5*(zp1(1,1,1)+zp1(1,1,2)))
      DO i=2,nx1-2
        DO j=2,ny1-2
          zpmin1=MIN(zpmin1,0.5*(zp1(i,j,1)+zp1(i,j,2)))
        END DO
      END DO

      snddelz = (zpmax-zpmin)/(lvlprof-2)
      DO k=1,lvlprof
        zsnd(k)= zpmin+(k-2)*snddelz
      END DO

      IF (bglopt == 4) THEN

        DO k=2,lvlprof
          temz1d3(k)=0.0
          temz1d4(k)=0.0
          temz1d5(k)=0.0
          temz1d6(k)=0.0
          temz1d7(k)=0.0
        END DO
        pbartop = 0.0

        DO i=2,nx-2
          DO j=2,ny-2

            DO k=2,nz-1
              temz1d1(k)=ptbar(i,j,k)
              temz1d2(k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
            END DO
            CALL inte1d(temz1d1(2),temz1d2(2),nz-2,                     &
                      ptsnd(2),zsnd(2),lvlprof-1)

            DO k=2,nz-1
              temz1d1(k)=qvbar(i,j,k)
            END DO
            CALL inte1d(temz1d1(2),temz1d2(2),nz-2,                     &
                      qvsnd(2),zsnd(2),lvlprof-1)

            DO k=2,nz-1
              temz1d1(k)=(ubar(i,j,k)+ubar(i+1,j,k))*0.5
            END DO
            CALL inte1d(temz1d1(2),temz1d2(2),nz-2,                     &
                      ubsnd(2),zsnd(2),lvlprof-1)

            DO k=2,nz-1
              temz1d1(k)=(vbar(i,j,k)+vbar(i,j+1,k))*0.5
            END DO
            CALL inte1d(temz1d1(2),temz1d2(2),nz-2,                     &
                      vbsnd(2),zsnd(2),lvlprof-1)

            DO k=2,nz-1
              temz1d1(k)=ALOG(pbar(i,j,k))
            END DO
            CALL inte1d(temz1d1(2),temz1d2(2),nz-2,                     &
                      tem,zsnd(lvlprof-1),1)

            pbartop = pbartop+tem

            DO k=2,lvlprof
              IF(zsnd(k) >= temz1d2(2)-1.0E-5) THEN
                 ! Add only points above ground
                temz1d3(k)=temz1d3(k)+1
                temz1d4(k)=temz1d4(k)+ptsnd(k)
                temz1d5(k)=temz1d5(k)+qvsnd(k)
                temz1d6(k)=temz1d6(k)+ubsnd(k)
                temz1d7(k)=temz1d7(k)+vbsnd(k)
              END IF
            END DO

          END DO
        END DO
!
!  Obtain averages
!
        DO k=2,lvlprof
          ptsnd(k)=temz1d4(k)/temz1d3(k)
          qvsnd(k)=temz1d5(k)/temz1d3(k)
          ubsnd(k)=temz1d6(k)/temz1d3(k)
          vbsnd(k)=temz1d7(k)/temz1d3(k)
        END DO

        pbartop=EXP(pbartop/temz1d3(lvlprof-1))
!
!  Integrate hydrostatic equation down from top
!
        DO k=2,lvlprof
          temz1d1(k) = ptsnd(k)*(1.0+rvdrd*qvsnd(k))/(1.0+qvsnd(k))
        END DO

        psnd(lvlprof-1) = ( pbartop/p0 )**(rddcp)

        DO k=lvlprof-2,2,-1
          psnd(k)=psnd(k+1)-g*(zsnd(k)-zsnd(k+1))                       &
                       /(cp*0.5*(temz1d1(k)+temz1d1(k+1)))
        END DO
        DO k=lvlprof,lvlprof
          psnd(k)=psnd(k-1)-g*(zsnd(k)-zsnd(k-1))                       &
                       /(cp*0.5*(temz1d1(k)+temz1d1(k-1)))
        END DO

        DO k = 2,lvlprof
          psnd(k) = psnd(k)**(1/rddcp) * p0
        END DO

        IF( zpmin1 < zpmin ) THEN  ! Extend 1d sounding to zpmin1 using
                                   ! constant temperature lapse rate

          dtdz = -0.0068  ! Minus lapse rate
          zsnd(1)=zpmin1
          delz = zpmin1 - zpmin

          tbark1 = ptsnd(2) * (psnd(2)/p0)**rddcp
          tbark  = tbark1+dtdz*delz

          lnpbar = ALOG(psnd(2))-g*2.0*delz/(rd*(tbark1+tbark))
          psnd(1) = EXP(lnpbar)

          ptsnd(1)=tbark*(p0/psnd(1))**rddcp
          qvsnd(1)=qvsnd(2)
          ubsnd(1)=ubsnd(2)
          vbsnd(1)=vbsnd(2)

        ELSE  ! Set for safty, won't be used.

          psnd(1)=psnd(2)
          ptsnd(1)=ptsnd(2)
          qvsnd(1)=qvsnd(2)
          ubsnd(1)=ubsnd(2)
          vbsnd(1)=vbsnd(2)

        END IF

      ELSE ! bglopt=5
!
!  Use extmnsnd (from ext2arps) to get average sounding
!
        ! tem1 - zp at zone center
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              tem1(i,j,k) = 0.5*(zp(i,j,k)+zp(i,j,k+1))
            END DO
          END DO
        END DO
        CALL edgfill(tem1,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)

        ! tem2 = log(pbar)
        ! tem3 = t
        ! tem4 = rhstar
        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              tem2(i,j,k) = ALOG(pbar(i,j,k))                      ! log pressure
              tem3(i,j,k) = (ptprt(i,j,k)+ptbar(i,j,k))                 &
                         *(((pprt(i,j,k)+pbar(i,j,k))/p0)**rddcp)  ! temperature
              qvsat=f_qvsat( pprt(i,j,k)+pbar(i,j,k), tem3(i,j,k) )
              tem4(i,j,k)=SQRT(AMAX1(0.,(rhmax-(qv(i,j,k)/qvsat))))
            END DO
          END DO
        END DO

        DO j=1,ny-1
          DO i=1,nx-1
            xs_2d(i,j) = xs(i)
            ys_2d(i,j) = ys(j)
          END DO
        END DO

        i1mn = MAX(1,INT(xorig1/dx1-1.))
        i1mx = MIN(nx-1,i1mn+nx1+1)
        j1mn = MAX(1,INT(yorig1/dy1-1.))
        j1mx = MIN(ny-1,j1mn+ny1+1)

        CALL extmnsnd(nx1,ny1,nx,ny,nz,lvlprof,1,                       &
                        i1mn,i1mx,j1mn,j1mx,2,nz-1,                     &
                        xs1,ys1,xs_2d,ys_2d,                            &
                        tem1,tem2,tem3,tem4,                            &
                        u,v,   & ! note: not vital that u & v be staggered
                        zsnd,temz1d1,psnd,temz1d2,ptsnd,temz1d3,        &
                        qvsnd,temz1d4,ubsnd,vbsnd)

      END IF

      WRITE(6,'(/a/)')'    Reconstructed horizontal average sounding:'
      DO k = lvlprof,1,-1
        WRITE(6,'(a,i4,4f10.2,f10.6)') 'k,z,p,pt,t,qv=',                &
                k,zsnd(k),psnd(k),ptsnd(k),                             &
                ptsnd(k)*(psnd(k)/p0)**rddcp,qvsnd(k)
      END DO

    END IF
!
!-----------------------------------------------------------------------
!
!  Intepolate 3D arrays at u-point
!
!-----------------------------------------------------------------------
!
    CALL avgz(zp1, 0 , nx1,ny1,nz1,1,nx1-1,1,ny1-1,1,nz1-1, tem11 )
    CALL avgsu(tem11,nx1,ny1,nz1,1,ny1-1,1,nz1-1, tem21, tem11)

    CALL avgz(zp, 0 , nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, tem1 )
    CALL avgsu(tem1,nx,ny,nz,1,ny-1,1,nz-1, tem2, tem1)

    CALL intrpxy3d(tem2,nx,1,nx,ny,1,ny-1,nz,1,nz-1,                    &
                   wgtux,iux,wgtsy,jsy,intrphopt,                       &
                   temx1y1zb,nx1,1,nx1, ny1,1,ny1-1, temx1yz)

    DO i=1,nx1
      DO j=1,ny1-1
        k0 = nz-2
        DO k=nz1-1,1,-1
          z0 = tem21(i,j,k)
          DO kk=k0,1,-1
            IF(z0 >= temx1y1zb(i,j,kk)) THEN
              k0 = kk
              EXIT
            END IF
          END DO
!          1115      CONTINUE
          IF(intrpvopt == 1) kz(i,j,k)=MAX(1,MIN(k0,nz-2))
          IF(intrpvopt == 2) kz(i,j,k)=MAX(2,MIN(k0,nz-2))
        END DO
      END DO
    END DO

    IF( intrpvopt == 1) THEN
      DO i=1,nx1
        DO j=1,ny1-1
          DO k=1,nz1-1
            wgtz(i,j,k,1)=(temx1y1zb(i,j,kz(i,j,k)+1)-tem21(i,j,k))/    &
                     (temx1y1zb(i,j,kz(i,j,k)+1)-temx1y1zb(i,j,kz(i,j,k)))
          END DO
        END DO
      END DO
    ELSE
      DO i=1,nx1
        DO j=1,ny1-1
          DO k=1,nz1-1
            d=tem21(i,j,k)
            a=temx1y1zb(i,j,kz(i,j,k)-1)
            b=temx1y1zb(i,j,kz(i,j,k))
            c=temx1y1zb(i,j,kz(i,j,k)+1)
            wgtz(i,j,k,1)= (d-b)*(d-c)/((a-b)*(a-c))
            wgtz(i,j,k,2)= (d-a)*(d-c)/((b-a)*(b-c))
            wgtz(i,j,k,3)= (d-a)*(d-b)/((c-a)*(c-b))
          END DO
        END DO
      END DO
    END IF

!    IF( nfile == 1) THEN
    IF(nfile == 1 .or. noutgrds > 0) THEN

      IF( redo_base_state == 0 ) THEN ! get base-state from original arrays

        CALL intrpxyz3d(ubar,nx,1,nx,ny,1,ny-1,nz,1,nz-1,               &
                        wgtux,iux,wgtsy,jsy,wgtz,kz,                    &
                        intrphopt,intrpvopt,                            &
                        ubar1,nx1,1,nx1, ny1,1,ny1-1,nz1,1,nz1-1,       &
                        temx1yz,temx1y1z)
        DO i=1,nx1
          DO j=1,ny1-1
            DO k=1,nz1-1
              IF(tem21(i,j,k) < temx1y1zb(i,j,2)) THEN ! below input grid ground
                IF(bglopt == 2) THEN
                  ubar1(i,j,k)=temx1y1z(i,j,2)
                ELSE IF(bglopt == 3) THEN
                  ubar1(i,j,k)=0.0
                END IF
              END IF
            END DO
          END DO
        END DO

      ELSE IF( redo_base_state == 1 ) THEN ! reconstruct base-state from 1d soundings

        DO k=1,nz1-1
          DO j=1,ny1-1
            DO i=1,nx1
              kk = MAX(1,MIN(lvlprof-1,                                 &
                   INT((tem21(i,j,k)-zsnd(2))/snddelz)+2))
              IF(tem21(i,j,k) < zsnd(2)) kk = 1
              tem = (tem21(i,j,k)-zsnd(kk))/(zsnd(kk+1)-zsnd(kk))
              ubar1(i,j,k)= ubsnd(kk)+(ubsnd(kk+1)-ubsnd(kk))*tem
            END DO
          END DO
        END DO

      END IF

    END IF

    IF (ntagopt == 1) THEN ! adjust tem21 (heights) so that surface wind is
                           ! same as input surface wind (by matching tem21
                           ! to the input height, temx1y1zb, near the surface).
      DO i=1,nx1
        DO j=1,ny1-1
          hght0 = tem21(i,j,2)
          dhght = hght0 - temx1y1zb(i,j,2)
          IF (dhght > 0) THEN
            DO k=1,nz1-1
              fac = MAX(0., tem21(i,j,k) - hght0)/aghght
              IF (fac > 1) EXIT
              tem21(i,j,k) = fac*fac*tem21(i,j,k)  &
                             + (1.-fac*fac)*(tem21(i,j,k)-dhght)
            END DO
          ENDIF
        END DO
      END DO

      DO i=1,nx1
        DO j=1,ny1-1
          k0 = nz-2
          DO k=nz1-1,1,-1
            z0 = tem21(i,j,k)
            DO kk=k0,1,-1
              IF(z0 >= temx1y1zb(i,j,kk)) THEN
                k0 = kk
                EXIT
              END IF
            END DO
            IF(intrpvopt == 1) kz(i,j,k)=MAX(1,MIN(k0,nz-2))
            IF(intrpvopt == 2) kz(i,j,k)=MAX(2,MIN(k0,nz-2))
          END DO
        END DO
      END DO

      IF( intrpvopt == 1) THEN
        DO i=1,nx1
          DO j=1,ny1-1
            DO k=1,nz1-1
              wgtz(i,j,k,1)=(temx1y1zb(i,j,kz(i,j,k)+1)-tem21(i,j,k))/    &
                       (temx1y1zb(i,j,kz(i,j,k)+1)-temx1y1zb(i,j,kz(i,j,k)))
            END DO
          END DO
        END DO
      ELSE
        DO i=1,nx1
          DO j=1,ny1-1
            DO k=1,nz1-1
              d=tem21(i,j,k)
              a=temx1y1zb(i,j,kz(i,j,k)-1)
              b=temx1y1zb(i,j,kz(i,j,k))
              c=temx1y1zb(i,j,kz(i,j,k)+1)
              wgtz(i,j,k,1)= (d-b)*(d-c)/((a-b)*(a-c))
              wgtz(i,j,k,2)= (d-a)*(d-c)/((b-a)*(b-c))
              wgtz(i,j,k,3)= (d-a)*(d-b)/((c-a)*(c-b))
            END DO
          END DO
        END DO
      END IF

    END IF

    CALL intrpxyz3d(u,nx,1,nx,ny,1,ny-1,nz,1,nz-1,                      &
                    wgtux,iux,wgtsy,jsy,wgtz,kz,                        &
                    intrphopt,intrpvopt,                                &
                    u1,nx1,1,nx1, ny1,1,ny1-1,nz1,1,nz1-1,              &
                    temx1yz,temx1y1z)

    IF( npoint_below_ground /= 0 ) THEN ! Set values below input grid ground

      DO i=1,nx1
        DO j=1,ny1-1
          DO k=nz1-2,1,-1
            IF(tem21(i,j,k) < temx1y1zb(i,j,2)) THEN ! below input grid ground
              IF(bglopt == 2 .OR. bglopt == 4 .OR. bglopt == 5) THEN
                u1   (i,j,k)=temx1y1z(i,j,2)
              ELSE IF(bglopt == 3) THEN
                u1   (i,j,k)=misvalue
              END IF
            END IF
          END DO
        END DO
      END DO

    END IF
!
!-----------------------------------------------------------------------
!
!  Intepolate 3D arrays at v-point
!
!-----------------------------------------------------------------------
!
    CALL avgz(zp1, 0, nx1,ny1,nz1,1,nx1-1,1,ny1-1,1,nz1-1, tem11 )
    CALL avgsv(tem11,nx1,ny1,nz1,1,nx1-1,1,nz1-1, tem21, tem11)

    CALL avgz(zp, 0, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, tem1 )
    CALL avgsv(tem1,nx,ny,nz,1,nx-1,1,nz-1, tem2, tem1)

    CALL intrpxy3d(tem2,nx,1,nx-1,ny,1,ny,nz,1,nz-1,                    &
                   wgtsx,isx,wgtvy,jvy,intrphopt,                       &
                   temx1y1zb,nx1,1,nx1-1, ny1,1,ny1, temx1yz)

    DO i=1,nx1-1
      DO j=1,ny1

        k0=nz-2
        DO k=nz1-1,1,-1
          z0 = tem21(i,j,k)
          DO kk=k0,1,-1
            IF(z0 >= temx1y1zb(i,j,kk)) THEN
              k0 = kk
              EXIT
            END IF
          END DO
!          1215      CONTINUE
          IF(intrpvopt == 1) kz(i,j,k)=MAX(1,MIN(k0,nz-2))
          IF(intrpvopt == 2) kz(i,j,k)=MAX(2,MIN(k0,nz-2))
        END DO

      END DO
    END DO

    IF( intrpvopt == 1) THEN
      DO i=1,nx1-1
        DO j=1,ny1
          DO k=1,nz1-1
            wgtz(i,j,k,1)=(temx1y1zb(i,j,kz(i,j,k)+1)-tem21(i,j,k))/    &
                      (temx1y1zb(i,j,kz(i,j,k)+1)-temx1y1zb(i,j,kz(i,j,k)))
          END DO
        END DO
      END DO
    ELSE
      DO i=1,nx1-1
        DO j=1,ny1
          DO k=1,nz1-1
            d=tem21(i,j,k)
            a=temx1y1zb(i,j,kz(i,j,k)-1)
            b=temx1y1zb(i,j,kz(i,j,k))
            c=temx1y1zb(i,j,kz(i,j,k)+1)
            wgtz(i,j,k,1)= (d-b)*(d-c)/((a-b)*(a-c))
            wgtz(i,j,k,2)= (d-a)*(d-c)/((b-a)*(b-c))
            wgtz(i,j,k,3)= (d-a)*(d-b)/((c-a)*(c-b))
          END DO
        END DO
      END DO
    END IF

!    IF( nfile == 1) THEN
    IF(nfile == 1 .or. noutgrds > 0) THEN

      IF( redo_base_state == 0 ) THEN ! get base-state from original arrays

        CALL intrpxyz3d(vbar,nx,1,nx-1,ny,1,ny,nz,1,nz-1,               &
                        wgtsx,isx,wgtvy,jvy,wgtz,kz,                    &
                        intrphopt,intrpvopt,                            &
                        vbar1,nx1,1,nx1-1, ny1,1,ny1,nz1,1,nz1-1,       &
                        temx1yz,temx1y1z)
        DO i=1,nx1-1
          DO j=1,ny1
            DO k=nz1-2,1,-1
              IF(tem21(i,j,k) < temx1y1zb(i,j,2)) THEN ! below input grid ground
                IF(bglopt == 2) THEN
                  vbar1(i,j,k)=temx1y1z(i,j,2)
                ELSE IF(bglopt == 3) THEN
                  vbar1(i,j,k)=0.0
                END IF
              END IF
            END DO
          END DO
        END DO

      ELSE IF( redo_base_state == 1 ) THEN ! reconstruct base-state from 1d soundings

        DO k=1,nz1-1
          DO i=1,nx1-1
            DO j=1,ny1
              kk = MAX(1,MIN(lvlprof-1,                                 &
                   INT((tem21(i,j,k)-zsnd(2))/snddelz)+2))
              IF(tem21(i,j,k) < zsnd(2)) kk = 1
              tem = (tem21(i,j,k)-zsnd(kk))/(zsnd(kk+1)-zsnd(kk))
              vbar1(i,j,k)= vbsnd(kk)+(vbsnd(kk+1)-vbsnd(kk))*tem
            END DO
          END DO
        END DO

      END IF

    END IF

    IF (ntagopt == 1) THEN ! adjust tem21 (heights) so that surface wind is
                           ! same as input surface wind (by matching tem21
                           ! to the input height, temx1y1zb, near the surface).
      DO i=1,nx1-1
        DO j=1,ny1
          hght0 = tem21(i,j,2)
          dhght = hght0 - temx1y1zb(i,j,2)
          IF (dhght > 0) THEN
            DO k=1,nz1-1
              fac = MAX(0., tem21(i,j,k) - hght0)/aghght
              IF (fac > 1) EXIT
              tem21(i,j,k) = fac*fac*tem21(i,j,k)  &
                             + (1.-fac*fac)*(tem21(i,j,k)-dhght)
            END DO
          ENDIF
        END DO
      END DO

      DO i=1,nx1-1
        DO j=1,ny1

          k0=nz-2
          DO k=nz1-1,1,-1
            z0 = tem21(i,j,k)
            DO kk=k0,1,-1
              IF(z0 >= temx1y1zb(i,j,kk)) THEN
                k0 = kk
                EXIT
              END IF
            END DO
            IF(intrpvopt == 1) kz(i,j,k)=MAX(1,MIN(k0,nz-2))
            IF(intrpvopt == 2) kz(i,j,k)=MAX(2,MIN(k0,nz-2))
          END DO

        END DO
      END DO

      IF( intrpvopt == 1) THEN
        DO i=1,nx1-1
          DO j=1,ny1
            DO k=1,nz1-1
              wgtz(i,j,k,1)=(temx1y1zb(i,j,kz(i,j,k)+1)-tem21(i,j,k))/    &
                        (temx1y1zb(i,j,kz(i,j,k)+1)-temx1y1zb(i,j,kz(i,j,k)))
            END DO
          END DO
        END DO
      ELSE
        DO i=1,nx1-1
          DO j=1,ny1
            DO k=1,nz1-1
              d=tem21(i,j,k)
              a=temx1y1zb(i,j,kz(i,j,k)-1)
              b=temx1y1zb(i,j,kz(i,j,k))
              c=temx1y1zb(i,j,kz(i,j,k)+1)
              wgtz(i,j,k,1)= (d-b)*(d-c)/((a-b)*(a-c))
              wgtz(i,j,k,2)= (d-a)*(d-c)/((b-a)*(b-c))
              wgtz(i,j,k,3)= (d-a)*(d-b)/((c-a)*(c-b))
            END DO
          END DO
        END DO
      END IF

    END IF

    CALL intrpxyz3d(v,nx,1,nx-1,ny,1,ny,nz,1,nz-1,                      &
                    wgtsx,isx,wgtvy,jvy,wgtz,kz,                        &
                    intrphopt,intrpvopt,                                &
                    v1,nx1,1,nx1-1, ny1,1,ny1,nz1,1,nz1-1,              &
                    temx1yz,temx1y1z)

    IF( npoint_below_ground /= 0 ) THEN ! Set values below input grid ground

      DO i=1,nx1-1
        DO j=1,ny1
          DO k=nz1-2,1,-1
            IF(tem21(i,j,k) < temx1y1zb(i,j,2)) THEN ! below input grid ground
              IF(bglopt == 2 .OR. bglopt == 4 .OR. bglopt == 5) THEN
                v1   (i,j,k)=temx1y1z(i,j,2)
              ELSE IF(bglopt == 3) THEN
                v1   (i,j,k)=misvalue
              END IF
            END IF
          END DO
        END DO
      END DO

    END IF
!
!-----------------------------------------------------------------------
!
!  Intepolate 3D arrays at w-point
!
!-----------------------------------------------------------------------
!
    CALL intrpxy3d(zp,nx,1,nx-1,ny,1,ny-1,nz,1,nz,                      &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   temx1y1zb,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)

    ! use tem21 instead of zp1 below
    tem21 = zp1

    IF (ntagopt == 1) THEN ! adjust tem21 (heights) so that surface wind is
                           ! same as input surface wind (by matching tem21
                           ! to the input height, temx1y1zb, near the surface).
      DO i=1,nx1-1
        DO j=1,ny1-1
          hght0 = tem21(i,j,2)
          dhght = hght0 - temx1y1zb(i,j,2)
          IF (dhght > 0) THEN
            DO k=1,nz1-1
              fac = MAX(0., tem21(i,j,k) - hght0)/aghght
              IF (fac > 1) EXIT
              tem21(i,j,k) = fac*fac*tem21(i,j,k)  &
                             + (1.-fac*fac)*(tem21(i,j,k)-dhght)
            END DO
          ENDIF
        END DO
      END DO

    END IF

    DO i=1,nx1-1
      DO j=1,ny1-1
        k0 = nz-1
        DO k=nz1,1,-1
          z0 = tem21(i,j,k)
          DO kk=k0,1,-1
            IF(z0 >= temx1y1zb(i,j,kk)) THEN
              k0 = kk
              EXIT
            END IF
          END DO
!          1315      CONTINUE

          IF(intrpvopt == 1) kz(i,j,k)=MAX(2,MIN(k0,nz-1))
          IF(intrpvopt == 2) kz(i,j,k)=MAX(3,MIN(k0,nz-1))
        END DO

      END DO
    END DO

    IF( intrpvopt == 1) THEN
      DO i=1,nx1-1
        DO j=1,ny1-1
          DO k=1,nz1
            wgtz(i,j,k,1)=(temx1y1zb(i,j,kz(i,j,k)+1)-tem21(i,j,k))/      &
                      (temx1y1zb(i,j,kz(i,j,k)+1)-temx1y1zb(i,j,kz(i,j,k)))
          END DO
        END DO
      END DO
    ELSE
      DO i=1,nx1-1
        DO j=1,ny1-1
          DO k=1,nz1
            d=tem21(i,j,k)
            a=temx1y1zb(i,j,kz(i,j,k)-1)
            b=temx1y1zb(i,j,kz(i,j,k))
            c=temx1y1zb(i,j,kz(i,j,k)+1)
            wgtz(i,j,k,1)= (d-b)*(d-c)/((a-b)*(a-c))
            wgtz(i,j,k,2)= (d-a)*(d-c)/((b-a)*(b-c))
            wgtz(i,j,k,3)= (d-a)*(d-b)/((c-a)*(c-b))
          END DO
        END DO
      END DO
    END IF

    CALL intrpxyz3d(w,nx,1,nx-1,ny,1,ny-1,nz,1,nz,                      &
                    wgtsx,isx,wgtsy,jsy,wgtz,kz,                        &
                    intrphopt,intrpvopt,                                &
                    w1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1,              &
                    temx1yz,temx1y1z)

    IF( npoint_below_ground /= 0 ) THEN ! Set values below input grid ground

      DO i=1,nx1-1
        DO j=1,ny1-1
          DO k=nz1-1,1,-1
            IF(zp1(i,j,k) < temx1y1zb(i,j,2)) THEN ! Blow input grid ground
              IF(bglopt == 2 .OR. bglopt == 4 .OR. bglopt == 5) THEN
                w1(i,j,k)=temx1y1z(i,j,2)
              ELSE IF(bglopt == 3) THEN
                w1(i,j,k)=misvalue
              END IF
            END IF
          END DO
        END DO
      END DO

    END IF

!-----------------------------------------------------------------------
!
!  Interpolate 3D arrays at scalar point
!
!-----------------------------------------------------------------------

    CALL avgz(zp1,0,nx1,ny1,nz1,1,nx1-1,1,ny1-1,1,nz1-1,tem21)

    CALL avgz(zp,0,nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem2)

    CALL intrpxy3d(tem2,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,                  &
                   wgtsx,isx,wgtsy,jsy,intrphopt,                       &
                   temx1y1zb,nx1,1,nx1-1, ny1,1,ny1-1, temx1yz)


    DO i=1,nx1-1
      DO j=1,ny1-1
        htrnx1y1(i,j)=temx1y1zb(i,j,2) ! Old terrain on new grid at scalar level
                                       ! e.g., z of first grid level AGL.
        DO k=1,nz1-1
          IF(tem21(i,j,k) >= htrnx1y1(i,j)) GO TO 1350
        END DO
        1350    CONTINUE
        ktrnx1y1(i,j)= k-1 ! The first scale-level below the first scale level
                           ! above input grid ground
      END DO
    END DO

    DO i=1,nx1-1
      DO j=1,ny1-1

        k0=nz-2
        DO k=nz1-1,1,-1
          z0 = tem21(i,j,k)
          DO kk=k0,1,-1
            IF(z0 >= temx1y1zb(i,j,kk) ) THEN
              k0 = kk
              EXIT
            END IF
          END DO
!          1415      CONTINUE
          IF(intrpvopt == 1) kz(i,j,k)=MAX(1,MIN(k0,nz-2))
          IF(intrpvopt == 2) kz(i,j,k)=MAX(2,MIN(k0,nz-2))
        END DO

      END DO
    END DO

    IF( intrpvopt == 1) THEN
      DO i=1,nx1-1
        DO j=1,ny1-1
          DO k=1,nz1-1
            wgtz(i,j,k,1)=(temx1y1zb(i,j,kz(i,j,k)+1)-tem21(i,j,k))/    &
                      (temx1y1zb(i,j,kz(i,j,k)+1)-temx1y1zb(i,j,kz(i,j,k)))
          END DO
        END DO
      END DO
    ELSE
      DO i=1,nx1-1
        DO j=1,ny1-1
          DO k=1,nz1-1
            d=tem21(i,j,k)
            a=temx1y1zb(i,j,kz(i,j,k)-1)
            b=temx1y1zb(i,j,kz(i,j,k))
            c=temx1y1zb(i,j,kz(i,j,k)+1)
            wgtz(i,j,k,1)= (d-b)*(d-c)/((a-b)*(a-c))
            wgtz(i,j,k,2)= (d-a)*(d-c)/((b-a)*(b-c))
            wgtz(i,j,k,3)= (d-a)*(d-b)/((c-a)*(c-b))
          END DO
        END DO
      END DO
    END IF

    IF(nfile == 1 .or. noutgrds > 0) THEN

      IF( redo_base_state == 0) THEN ! Get base-state from original arrays

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              tem1(i,j,k)=ALOG(pbar(i,j,k))
            END DO
          END DO
        END DO

        CALL intrpxyz3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,             &
                        wgtsx,isx,wgtsy,jsy,wgtz,kz,                    &
                        intrphopt,intrpvopt,                            &
                        pbar1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,     &
                        temx1yz,temx1y1z)

        DO k=1,nz1-1
          DO j=1,ny1-1
            DO i=1,nx1-1
              pbar1(i,j,k)=EXP(pbar1(i,j,k))
            END DO
          END DO
        END DO

        DO i=1,nx1-1
          DO j=1,ny1-1
            pbsfc(i,j)= EXP(temx1y1z(i,j,2))
          END DO
        END DO

        CALL intrpxyz3d(ptbar,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,            &
                        wgtsx,isx,wgtsy,jsy,wgtz,kz,                    &
                        intrphopt,intrpvopt,                            &
                        ptbar1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,    &
                        temx1yz,temx1y1z)

        DO i=1,nx1-1
          DO j=1,ny1-1
            ptbsfc(i,j)= temx1y1z(i,j,2)
          END DO
        END DO


        CALL intrpxyz3d(qvbar,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,            &
                        wgtsx,isx,wgtsy,jsy,wgtz,kz,                    &
                        intrphopt,intrpvopt,                            &
                        qvbar1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,    &
                        temx1yz,temx1y1z)

        DO i=1,nx1-1
          DO j=1,ny1-1
            qvbsfc(i,j)= temx1y1z(i,j,2)
          END DO
        END DO

        CALL intrpxyz3d(rhobar,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,           &
                        wgtsx,isx,wgtsy,jsy,wgtz,kz,                    &
                        intrphopt,intrpvopt,                            &
                        rhobar1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,   &
                        temx1yz,temx1y1z)

      ELSE IF( redo_base_state == 1 ) THEN ! reconstruct base-state from 1d soundings

        DO k=1,nz1-1
          DO i=1,nx1-1
            DO j=1,ny1-1

              kk = MAX(1,MIN(lvlprof-1,                                 &
                   INT((tem21(i,j,k)-zsnd(2))/snddelz)+2))
              IF(tem21(i,j,k) < zsnd(2)) kk = 1

              tem = (tem21(i,j,k)-zsnd(kk))/(zsnd(kk+1)-zsnd(kk))
              ptbar1(i,j,k)= ptsnd(kk)+(ptsnd(kk+1)-ptsnd(kk))*tem
              qvbar1(i,j,k)= qvsnd(kk)+(qvsnd(kk+1)-qvsnd(kk))*tem
              tema = ALOG(psnd (kk) )
              pbar1 (i,j,k)=EXP(tema+(ALOG(psnd(kk+1))-tema)*tem)
              rhobar1(i,j,k)=pbar1(i,j,k)/                              &
                            (rd*ptbar1(i,j,k)*(pbar1(i,j,k)/p0)**rddcp)

            END DO
          END DO
        END DO

        DO i=1,nx1-1
          DO j=1,ny1-1
            kk= MAX(1,MIN(lvlprof-1,                                    &
                INT((htrnx1y1(i,j)-zsnd(2))/snddelz)+2))
            IF(htrnx1y1(i,j) < zsnd(2)) k = 1
            tem = (htrnx1y1(i,j)-zsnd(kk))/(zsnd(kk+1)-zsnd(kk))
            ptbsfc(i,j)= ptsnd(kk)+(ptsnd(kk+1)-ptsnd(kk))*tem
            qvbsfc(i,j)= qvsnd(kk)+(qvsnd(kk+1)-qvsnd(kk))*tem
            pbsfc (i,j)= EXP(ALOG(psnd (kk))+                           &
                         (ALOG(psnd (kk+1))-ALOG(psnd (kk)))*tem)
          END DO
        END DO

      END IF

    END IF

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem2(i,j,k)=ALOG(pbar(i,j,k)+pprt(i,j,k))
        END DO
      END DO
    END DO

    CALL intrpxyz3d(tem2,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,                 &
                    wgtsx,isx,wgtsy,jsy,wgtz,kz,                        &
                    intrphopt,intrpvopt,                                &
                    pprt1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,         &
                    temx1yz,temx1y1z)

    DO k=1,nz1-1
      DO j=1,ny1-1
        DO i=1,nx1-1
          pprt1(i,j,k)=EXP(pprt1(i,j,k))-pbar1(i,j,k)
        END DO
      END DO
    END DO
    DO j=1,ny1-1
      DO i=1,nx1-1
        ppsfc(i,j)=EXP(temx1y1z(i,j,2))-pbsfc(i,j)
      END DO
    END DO

    IF ( redo_base_state == 1 ) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem2(i,j,k)=ptprt(i,j,k)+ptbar(i,j,k)
          END DO
        END DO
      END DO

      CALL intrpxyz3d(tem2,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,               &
                      wgtsx,isx,wgtsy,jsy,wgtz,kz,                      &
                      intrphopt,intrpvopt,                              &
                      ptprt1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,      &
                      temx1yz,temx1y1z)

      DO j=1,ny1-1
        DO i=1,nx1-1
          ptpsfc(i,j)=temx1y1z(i,j,2)-ptbsfc(i,j)
        END DO
      END DO

      DO k=1,nz1-1
        DO j=1,ny1-1
          DO i=1,nx1-1
            ptprt1(i,j,k)=ptprt1(i,j,k)-ptbar1(i,j,k)
          END DO
        END DO
      END DO

    ELSE

      CALL intrpxyz3d(ptprt,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,              &
                     wgtsx,isx,wgtsy,jsy,wgtz,kz,                       &
                     intrphopt,intrpvopt,                               &
                     ptprt1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,       &
                     temx1yz,temx1y1z)

      DO j=1,ny1-1
        DO i=1,nx1-1
          ptpsfc(i,j)=temx1y1z(i,j,2)
        END DO
      END DO

    END IF

    CALL intrpxyz3d(qv,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,                   &
                    wgtsx,isx,wgtsy,jsy,wgtz,kz,                        &
                    intrphopt,intrpvopt,                                &
                    qv1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,           &
                    temx1yz,temx1y1z)
    DO j=1,ny1-1
      DO i=1,nx1-1
        qvsfc(i,j)=temx1y1z(i,j,2)
      END DO
    END DO

    DO nq=1,nscalar

      CALL intrpxyz3d(qscalar(1,1,1,nq),nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,  &
                    wgtsx,isx,wgtsy,jsy,wgtz,kz,                        &
                    intrphopt,intrpvopt,                                &
                    qscalar1(1,1,1,nq),nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,  &
                    temx1yz,temx1y1z)
      DO j=1,ny1-1
        DO i=1,nx1-1
          qscalarsfc(i,j,nq)=temx1y1z(i,j,2)
        END DO
      END DO

    END DO

    CALL intrpxyz3d(tke,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,                  &
                    wgtsx,isx,wgtsy,jsy,wgtz,kz,                        &
                    intrphopt,intrpvopt,                                &
                    tke1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,          &
                    temx1yz,temx1y1z)
    DO j=1,ny1-1
      DO i=1,nx1-1
        tkesfc(i,j)=temx1y1z(i,j,2)
      END DO
    END DO

    CALL intrpxyz3d(kmh,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,                  &
                    wgtsx,isx,wgtsy,jsy,wgtz,kz,                        &
                    intrphopt,intrpvopt,                                &
                    kmh1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,          &
                    temx1yz,temx1y1z)
    DO j=1,ny1-1
      DO i=1,nx1-1
        kmhsfc(i,j)=temx1y1z(i,j,2)
      END DO
    END DO

    CALL intrpxyz3d(kmv,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,                  &
                    wgtsx,isx,wgtsy,jsy,wgtz,kz,                        &
                    intrphopt,intrpvopt,                                &
                    kmv1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,          &
                    temx1yz,temx1y1z)
    DO j=1,ny1-1
      DO i=1,nx1-1
        kmvsfc(i,j)=temx1y1z(i,j,2)
      END DO
    END DO

    CALL intrpxyz3d(radfrc,nx,1,nx-1,ny,1,ny-1,nz,1,nz-1,               &
                    wgtsx,isx,wgtsy,jsy,wgtz,kz,                        &
                    intrphopt,intrpvopt,                                &
                    radfrc1,nx1,1,nx1-1, ny1,1,ny1-1,nz1,1,nz1-1,       &
                    temx1yz,temx1y1z)

    DO j=1,ny1-1
      DO i=1,nx1-1
        radsfc(i,j)=temx1y1z(i,j,2)
      END DO
    END DO
!
!  Reset extrapolated values below ground for certain options
!
    IF( npoint_below_ground /= 0 ) THEN ! Set values below input grid ground

      IF( bglopt == 2 ) THEN ! Constant extension below ground

        DO i=1,nx1-1
          DO j=1,ny1-1
            DO k=1,ktrnx1y1(i,j) ! only below ground
              ptprt1(i,j,k) = ptpsfc(i,j)
              pprt1 (i,j,k) = ppsfc (i,j)
              qv1   (i,j,k) = qvsfc (i,j)
              DO nq=1,nscalar
                qscalar1(i,j,k,nq) = qscalarsfc(i,j,nq)
              END DO
              tke1  (i,j,k) = tkesfc(i,j)
              kmh1  (i,j,k) = kmhsfc(i,j)
              kmv1  (i,j,k) = kmhsfc(i,j)

              qvbar1(i,j,k) = qvbsfc(i,j)
              ptbar1(i,j,k) = ptbsfc(i,j)
              pbar1 (i,j,k) = pbsfc (i,j)
              rhobar1(i,j,k)= pbsfc(i,j)/                               &
                             (rd*ptbsfc(i,j)*(pbsfc(i,j)/p0)**rddcp)
              radfrc1(i,j,k)= radsfc(i,j)
            END DO
          END DO
        END DO

      ELSE IF( bglopt == 3 ) THEN ! Set to a missing value

        DO i=1,nx1-1
          DO j=1,ny1-1
            DO k=1,ktrnx1y1(i,j) ! only below ground
              ptprt1(i,j,k) = misvalue
              pprt1 (i,j,k) = misvalue
              qvbar1(i,j,k) = 0.0
              qv1   (i,j,k) = misvalue
              DO nq=1,nscalar
                qscalar1(i,j,k,nq) = misvalue
              END DO
              tke1  (i,j,k) = misvalue
              kmh1  (i,j,k) = misvalue
              kmv1  (i,j,k) = misvalue
              ptbar1(i,j,k) = 0.0
              pbar1 (i,j,k) = 0.0
              rhobar1(i,j,k)= misvalue
              radfrc1(i,j,k)= misvalue
            END DO
          END DO
        END DO
!
      ELSE IF( bglopt == 4  .OR. bglopt == 5) THEN
         ! Using a constant lapse rate for temperature
         ! and then the hydrostatic balance for pressure

        DO i=1,nx1-1
          DO j=1,ny1-1
            DO k=1,ktrnx1y1(i,j) ! only below ground
              qv1   (i,j,k) = qvsfc (i,j)
              DO nq=1,nscalar
                qscalar1(i,j,k,nq) = qscalarsfc(i,j,nq)
              END DO
              tke1  (i,j,k) = tkesfc(i,j)
              kmh1  (i,j,k) = kmhsfc(i,j)
              kmv1  (i,j,k) = kmhsfc(i,j)
              radfrc1(i,j,k)= radsfc(i,j)

              qvbar1(i,j,k) = qvbsfc(i,j)

            END DO
          END DO
        END DO

        IF( redo_base_state == 1 ) THEN ! There are points below ground

          dtdz = -0.0068  ! Minus lapse rate

          DO i=1,nx1-1
            DO j=1,ny1-1

              k=ktrnx1y1(i,j)  ! First level below input grid ground
              IF( k >= 1) THEN
                delz = tem21(i,j,k)-htrnx1y1(i,j)
                ttotk1 = (ptbsfc(i,j)+ptpsfc(i,j)) *                    &
                        ((pbsfc (i,j)+ppsfc (i,j))/p0)**rddcp
                ttotk  = ttotk1+dtdz*delz
                lnptot = ALOG(pbsfc (i,j)+ppsfc (i,j))                  &
                         -g*2.0*delz/(rd*(ttotk1+ttotk))
                pprt1(i,j,k) = EXP(lnptot)-pbar1(i,j,k)
                ptprt1(i,j,k)= ttotk*(p0/(pbar1(i,j,k)+pprt1(i,j,k)))   &
                               **rddcp-ptbar1(i,j,k)
              END IF

              DO k=ktrnx1y1(i,j)-1,1,-1 ! other below-ground levels
                delz = tem21(i,j,k)-tem21(i,j,k+1)
                ttotk1 = (ptbar1(i,j,k+1)+ptprt1(i,j,k+1)) *            &
                        ((pbar1 (i,j,k+1)+pprt1 (i,j,k+1))/p0)**rddcp
                ttotk  = ttotk1+dtdz*delz
                lnptot = ALOG(pbar1(i,j,k+1)+pprt1(i,j,k+1))            &
                         -g*2.0*delz/(rd*(ttotk1+ttotk))
                pprt1(i,j,k) = EXP(lnptot)-pbar1(i,j,k)
                ptprt1(i,j,k)= ttotk*(p0/(pbar1(i,j,k)+pprt1(i,j,k)))   &
                               **rddcp-ptbar1(i,j,k)
              END DO
            END DO
          END DO

        END IF

      END IF

    END IF ! Check on npoint_below_ground

    800   CONTINUE

    xgrdorg = xorig1
    ygrdorg = yorig1

    CALL xytoll(1,1,xgrdorg,ygrdorg,origlat,origlon)
    CALL setorig( 2, origlat, origlon) ! set up the model origin on new grid

    IF(nfile == 1 .or. noutgrds > 0) THEN

      DO i=1,nx1
        x1_out(i) = x1(i)-xgrdorg
      END DO

      DO j=1,ny1
        y1_out(j) = y1(j)-ygrdorg
      END DO

      WRITE (*,*) "ctrlat & ctrlon:",ctrlat1,ctrlon1

      WRITE(6,'(/2x,a)')'For the output grid:'

      WRITE(6,'(2x,a,f15.7)')                                               &
          'The latitude  of the output grid center, ctrlat1=',ctrlat1
      WRITE(6,'(2x,a,f15.7/)')                                              &
          'The longitude of the output grid center, ctrlon1=',ctrlon1
      WRITE(6,'(2x,a,f15.7)')                                               &
          'The x coordinate of the output grid center, xctr1=',xctr1
      WRITE(6,'(2x,a,f15.7)')                                               &
          'The y coordinate of the output grid center, yctr1=',yctr1

      IF( ABS(xctr1-xctr1_old) > 0.0001*dx .OR.                             &
            ABS(yctr1-yctr1_old) > 0.0001*dy ) THEN
        WRITE(6,'(/1x,a/)')                                                 &
            '################################################################'
        WRITE(6,'(/2x,a/2x,a)')                                             &
            'Note that xctr1 and/or yctr1 had been changed by the program', &
            'due to snap_to_grid option. Check earlier messages.'
        WRITE(6,'(/1x,a/)')                                                 &
            '################################################################'
      END IF

      WRITE(6,'(/2x,a/2x,a,2f15.7,a)')                                      &
          'The SW corner (i,j)=(2,2) of the grid is located at ',           &
          '(',xgrdorg,ygrdorg,') of the input grid.'

    END IF


!-----------------------------------------------------------------------
!
!  Copy grid information to the common block variables so that output
!  routines have the correct values for the new grid.
!
!-----------------------------------------------------------------------

    ! save global setting before setting them to new grid values
    ctrlat_sv = ctrlat
    ctrlon_sv = ctrlon
    latitud_sv = latitud
    dx_sv = dx
    dy_sv = dy
    dz_sv = dz
    dzmin_sv = dzmin
    strhopt_sv = strhopt
    zrefsfc_sv = zrefsfc
    dlayer1_sv = dlayer1
    dlayer2_sv = dlayer2
    zflat_sv = zflat
    strhtune_sv = strhtune
    nstyp_sv = nstyp

    ctrlat = ctrlat1
    ctrlon = ctrlon1
    latitud = ctrlat

    dx = dx1
    dy = dy1
    dz = dz1

    dzmin = dzmin1
    strhopt = strhopt1
    zrefsfc = zrefsfc1
    dlayer1 = dlayer11
    dlayer2 = dlayer21
    zflat = zflat1
    strhtune = strhtune1

    nstyp = nstyp1

    WRITE (cmnt(nocmnt),'(a,i4,a,i4,a,i4)')                             &
        ' nx =',nx1,', ny =',ny1,', nz =',nz1
!
!-----------------------------------------------------------------------
!
!  Print out the max/min of output varaibles.
!
!-----------------------------------------------------------------------
!
    WHERE(tke1 <= 0.0) tke1 = 0.0

    WRITE(6,'(/1x,a/)')                                                 &
        'Min. and max. of input data interpolated to the new grid:'

    CALL a3dmax0(x1_out,1,nx1,1,nx1,1,1,1,1,1,1,1,1, amax,amin)
    WRITE(6,'(/1x,2(a,e13.6))') 'xmin    = ', amin,',  xmax    =',amax

    CALL a3dmax0(y1_out,1,ny1,1,ny1,1,1,1,1,1,1,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'ymin    = ', amin,',  ymax    =',amax

    CALL a3dmax0(z1,1,nz1,1,nz1,1,1,1,1,1,1,1,1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'zmin    = ', amin,',  zmax    =',amax

    CALL a3dmax0(zp1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1,           &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'zpmin   = ', amin,', zpmax    =',amax

    CALL a3dmax0(ubar1,1,nx1,1,nx1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,         &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'ubarmin = ', amin,',  ubarmax =',amax

    CALL a3dmax0(vbar1,1,nx1,1,nx1-1,1,ny1,1,ny1,1,nz1,1,nz1-1,         &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'vbarmin = ', amin,',  vbarmax =',amax

    CALL a3dmax0(ptbar1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,      &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'ptbarmin= ', amin,',  ptbarmax=',amax

    CALL a3dmax0(pbar1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,       &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'pbarmin = ', amin,',  pbarmax =',amax

    CALL a3dmax0(rhobar1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,     &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'rhobarmin=', amin,', rhobarmax=',amax

    CALL a3dmax0(qvbar1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,      &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'qvbarmin= ', amin,',  qvbarmax=',amax

    CALL a3dmax0(u1,1,nx1,1,nx1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,            &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax

    CALL a3dmax0(v1,1,nx1,1,nx1-1,1,ny1,1,ny1,1,nz1,1,nz1-1,            &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax

    CALL a3dmax0(w1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1,            &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax

    CALL a3dmax0(ptprt1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,      &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'ptprtmin= ', amin,',  ptprtmax=',amax

    CALL a3dmax0(pprt1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,       &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'pprtmin = ', amin,',  pprtmax =',amax

    CALL a3dmax0(qv1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,         &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'qvmin   = ', amin,',  qvmax   =',amax

    DO nq=1,nscalar
      CALL a3dmax0(qscalar1(1,1,1,nq),1,nx1,1,nx1-1,1,ny1,1,ny1-1,         &
                   1,nz1,1,nz1-1,amax,amin)

       WRITE(6,'(1x,2(a,e13.6))') TRIM(qnames(nq))//'min   = ', amin,    &
                           ',  '//TRIM(qnames(nq))//'max   =',  amax
    END DO

    CALL a3dmax0(tke1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,        &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'tkemin  = ', amin,',  tkemax  =',amax

    CALL a3dmax0(kmh1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,        &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'kmhmin  = ', amin,',  kmhmax  =',amax

    CALL a3dmax0(kmv1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,        &
                amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'kmvmin  = ', amin,',  kmvmax  =',amax

    CALL a3dmax0(raing1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'raingmin= ', amin,',  raingmax=',amax

    CALL a3dmax0(rainc1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'raincmin= ', amin,',  raincmax=',amax

    CALL a3dmax0(prcrate1(1,1,1),1,nx1,1,nx1-1,1,ny1,1,ny1-1,           &
                 1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'prcr1min= ', amin,',  prcr1max=',amax

    CALL a3dmax0(prcrate1(1,1,2),1,nx1,1,nx1-1,1,ny1,1,ny1-1,           &
                 1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'prcr2min= ', amin,',  prcr2max=',amax

    CALL a3dmax0(prcrate1(1,1,3),1,nx1,1,nx1-1,1,ny1,1,ny1-1,           &
                 1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'prcr3min= ', amin,',  prcr3max=',amax

    CALL a3dmax0(prcrate1(1,1,4),1,nx1,1,nx1-1,1,ny1,1,ny1-1,           &
                 1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'prcr4min= ', amin,',  prcr4max=',amax

    DO iss = 0, nstyp1

      CALL a3dmax0(tsoil1(1,1,1,iss),1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,  &
                   amax,amin)
      WRITE(6,'(1x,2(a,e13.6),a,i3)')                                   &
          'tsoil_tsfcmin = ', amin,',  tsoil_tsfcmax =',amax,' for soil type=',iss

      CALL a3dmax0(tsoil1(1,1,2,iss),1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1, &
                   amax,amin)
      WRITE(6,'(1x,2(a,e13.6),a,i3)')                                   &
          'tsoil_dpmin= ', amin,',  tsoil_dpmax=',amax,' for soil type=',iss

      CALL a3dmax0(qsoil1(1,1,1,iss),1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1, &
                   amax,amin)
      WRITE(6,'(1x,2(a,e13.6),a,i3)')                                   &
          'qsoil_wetsmin = ', amin,',  qsoil_wetsmax =',amax,' for soil type=',iss

      CALL a3dmax0(qsoil1(1,1,2,iss),1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1, &
                   amax,amin)
      WRITE(6,'(1x,2(a,e13.6),a,i3)')                                   &
          'qsoil_wetdmin = ', amin,',  qsoil_wetdmax =',amax,' for soil type=',iss

      CALL a3dmax0(wetcanp1(1,1,iss),                                   &
                   1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,amax,amin)
      WRITE(6,'(1x,2(a,e13.6),a,i3)')                                   &
          'wetcmin = ', amin,',  wetcmax =',amax,' for soil type=',iss

    END DO

    CALL a3dmax0(roufns1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,           &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'roufnmin= ', amin,', roufnmax =',amax

    CALL a3dmax0(veg1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,              &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'vegmin  = ', amin,',   vegmax =',amax

    CALL a3dmax0(radfrc1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,     &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'radfnmin= ', amin,', radfnmax =',amax

    CALL a3dmax0(radsw1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,            &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'radswmin= ', amin,', radswmax =',amax

    CALL a3dmax0(rnflx1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,            &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'rnflxmin= ', amin,', rnflxmax =',amax

    CALL a3dmax0(radswnet1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,         &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'radswnetmin= ', amin,', radswnetmax =',amax

    CALL a3dmax0(radlwin1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,          &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'radlwinmin= ', amin,', radlwinmax =',amax

    CALL a3dmax0(usflx1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,            &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'usflxmin= ', amin,', usflxmax =',amax

    CALL a3dmax0(vsflx1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,            &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'vsflxmin= ', amin,', vsflxmax =',amax

    CALL a3dmax0(ptsflx1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,           &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'ptflxmin= ', amin,', ptflxmax =',amax

    CALL a3dmax0(qvsflx1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,1,1,1,           &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'qvflxmin= ', amin,', qvflxmax =',amax

    tem11= 0.0        ! To be put in place of wbar
!
!-----------------------------------------------------------------------
!
!  Data dump of the model grid and base state arrays:
!
!  First find a unique name basdmpfn(1:lbasdmpf) for the grid and
!  base state array dump file
!
!  If grid/base state data has been written out once, skip
!  the following writing block. Also no need to write out
!  separate data for Savi3D dump. The same for GrADS dump.
!
!-----------------------------------------------------------------------
!
    runname = new_runname
    IF (noutgrds > 0) runname = name_grd(ng)

    CALL gtlfnkey(runname, lfnkey)

    IF(houtfmt /= 9 ) THEN

      IF( nfile > 1 ) GO TO 300 ! If done already, skip this part.

      CALL gtbasfn(runname(1:lfnkey),dirname,ldirnam,hdmpfmt,           &
                   1,0,basdmpfn,lbasdmpf)

      PRINT*
      PRINT*,'Output grid/base state file is ', basdmpfn(1:lbasdmpf)

      grdbas = 1      ! Dump out grd and base state arrays only

      CALL setcornerll( nx1,ny1,x1_out,y1_out )

      CALL dtadump(nx1,ny1,nz1,nzsoil1,nstyp1,hdmpfmt,nchout,           &
                   basdmpfn(1:lbasdmpf),grdbas,filcmprs,                &
                   u1,v1,w1,ptprt1,pprt1,qv1,qscalar1,                  &
                   tke1,kmh1,kmv1,                                      &
                   ubar1,vbar1,tem11,ptbar1,pbar1,rhobar1,qvbar1,       &
                   x1_out,y1_out,z1,zp1,zpsoil1,                        &
                   soiltyp1,stypfrct1,vegtyp1,lai1,roufns1,veg1,        &
                   tsoil1,qsoil1,wetcanp1,snowdpth1,                    &
                   raing1,rainc1,prcrate1,                              &
                   radfrc1,radsw1,rnflx1,radswnet1,radlwin1,            &
                   usflx1,vsflx1,ptsflx1,qvsflx1,                       &
                   tem21,wgtz,kz)


      300     CONTINUE

    END IF
!
!-----------------------------------------------------------------------
!
!  Then the time dependent fields:
!
!-----------------------------------------------------------------------
!
    IF( .NOT. (houtfmt == 9 .AND. vroutcnt == 1) ) THEN
!
!-----------------------------------------------------------------------
!
!  Reconstruct the file name using the specified directory name
!
!-----------------------------------------------------------------------
!
      CALL gtdmpfn(runname(1:lfnkey),dirname,                           &
                   ldirnam,curtim,hdmpfmt,1,0, hdmpfn, ldmpf)

    END IF

    IF ( exbchdfcompr > 4 ) rayklow = -1

    WRITE(6,'(a,a)') 'Writing t-dependent variable history dump ',      &
                      hdmpfn(1:ldmpf)
    grdbas = 0

    CALL dtadump(nx1,ny1,nz1,nzsoil1,nstyp1,hdmpfmt,nchout,             &
                 hdmpfn(1:ldmpf),grdbas,filcmprs,                       &
                 u1,v1,w1,ptprt1,pprt1,qv1,qscalar1,                    &
                 tke1,kmh1,kmv1,                                        &
                 ubar1,vbar1,tem11,ptbar1,pbar1,rhobar1,qvbar1,         &
                 x1_out,y1_out,z1,zp1,zpsoil1,                          &
                 soiltyp1,stypfrct1,vegtyp1,lai1,roufns1,veg1,          &
                 tsoil1,qsoil1,wetcanp1,snowdpth1,                      &
                 raing1,rainc1,prcrate1,                                &
                 radfrc1,radsw1,rnflx1,radswnet1,radlwin1,              &
                 usflx1,vsflx1,ptsflx1,qvsflx1,                         &
                 tem21,wgtz,kz)
!
!-----------------------------------------------------------------------
!
!  Write out soil model variable file
!
!-----------------------------------------------------------------------
!
    IF (sfcin == 1 .AND. soildmp > 0) THEN

      CALL cvttsnd( curtim, timsnd, tmstrln )

      soiloutfl = runname(1:lfnkey)//".soilvar."//timsnd(1:tmstrln)
      lfn = lfnkey + 9 + tmstrln

      IF( dirname /= ' ' ) THEN
        temchar = soiloutfl
        soiloutfl = dirname(1:ldirnam)//'/'//temchar
        lfn  = lfn + ldirnam + 1
      END IF

      !CALL fnversn(soiloutfl, lfn)

      !PRINT *, 'Write soil initial data to ',soiloutfl(1:lfn)

      CALL wrtsoil(nx1,ny1,nzsoil1,nstyp1, soiloutfl(1:lfn),            &
                   dx1,dy1,zpsoil1,                                     &
                   mapproj,trulat1,trulat2,trulon,sclfct,               &
                   ctrlat1,ctrlon1,                                     &
                   1,1,1,1,1,                                           &
                   tsoil1,qsoil1,wetcanp1,snowdpth1,soiltyp1)

      IF (soildmp == 1) CALL soilcntl(nx1,ny1,nzsoil1,zpsoil1,          &
                                      soiloutfl(1:lfn),                 &
                                      1,1,1,1,1, x1_out,y1_out)

    END IF       ! sfcin.eq.1

    IF(nfile == 1) THEN  ! Need to do this only once
!
!-----------------------------------------------------------------------
!
!  Write out terrain data
!
!-----------------------------------------------------------------------
!
      IF (terndmp > 0) THEN
        CALL getunit( nunit )

        ternfn = runname(1:lfnkey)//".trndata"
        lternfn = lfnkey + 8

        IF( dirname /= ' ' ) THEN

          temchar = ternfn
          ternfn = dirname(1:ldirnam)//'/'//temchar
          lternfn  = lternfn + ldirnam + 1

        END IF

        CALL fnversn(ternfn, lternfn )

        PRINT *, 'Write terrain data to ',ternfn(1:lternfn)

        CALL writtrn(nx1,ny1,ternfn(1:lternfn), dx1,dy1,                &
                   mapproj,trulat1,trulat2,trulon,sclfct,               &
                   ctrlat1,ctrlon1,hterain1)

        IF (terndmp == 1)                                               &
            CALL trncntl(nx1,ny1,ternfn(1:lternfn), x1_out,y1_out)

      END IF

!-----------------------------------------------------------------------
!
!  Write out surface property data file: sfcoutfl .
!
!-----------------------------------------------------------------------
!
      IF (landin == 1 .AND. sfcdmp > 0) THEN

        sfcoutfl = runname(1:lfnkey)//".sfcdata"
        lfn = lfnkey + 8

        IF( dirname /= ' ' ) THEN
          temchar = sfcoutfl
          sfcoutfl = dirname(1:ldirnam)//'/'//temchar
          lfn  = lfn + ldirnam + 1
        END IF

        CALL fnversn(sfcoutfl, lfn)

        PRINT *, 'Write surface property data in ',sfcoutfl(1:lfn)

        CALL wrtsfcdt(nx1,ny1,nstyp1,sfcoutfl(1:lfn), dx1,dy1,          &
                    mapproj,trulat1,trulat2,trulon,sclfct,              &
                    ctrlat1,ctrlon1,                                    &
                    1,1,1,1,1,0,                                        &
                    soiltyp1,stypfrct1,vegtyp1,lai1,roufns1,veg1,veg1)

        IF (sfcdmp == 1) CALL sfccntl(nx1,ny1, sfcoutfl(1:lfn),         &
                     1,1,1,1,1,0, x1_out,y1_out,tem11(1,1,1),           &
                     tem11(1,1,2))

      END IF       ! landin.eq.1

    END IF

    ctrlat = ctrlat_sv
    ctrlon = ctrlon_sv
    latitud = latitud_sv
    dx = dx_sv
    dy = dy_sv
    dz = dz_sv
    dzmin = dzmin_sv
    strhopt = strhopt_sv
    zrefsfc = zrefsfc_sv
    dlayer1 = dlayer1_sv
    dlayer2 = dlayer2_sv
    zflat = zflat_sv
    strhtune = strhtune_sv
    nstyp = nstyp_sv


    600 CONTINUE

    IF (noutgrds > 0 .and. ng < noutgrds) GOTO 500

    vroutcnt = 1 ! Variables have been written out at least once

    IF (hinfmt == 9) GO TO 102

  END DO

  9001  CONTINUE


  WRITE (6,'(/a,g13.3,a)') "Current memory allocation=",                &
               current_memory_use,' (Mw).'
  WRITE (6,'(a,g13.3,a/)') "Maxumum memory allocation=",                &
               max_memory_use,' (Mw).'

  gcpu_time = f_cputime() - cpu0
  WRITE (6,'(a,g13.3,a/)') 'Total CPU time used =', gcpu_time,' (s).'

  WRITE(6,'(a)') ' ==== Normal successful completion of ARPSTINTRP. ===='
  STOP

  9002  CONTINUE
  WRITE(6,'(1x,a,i2,/1x,a)')                                            &
      'Data read was unsuccessful. ireturn =', ireturn,                 &
      'Job stopped in ARPSINTRP.'

  STOP 9002
END PROGRAM arpsintrp
