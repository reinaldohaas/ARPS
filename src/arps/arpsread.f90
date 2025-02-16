PROGRAM arpsread
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Sample program to read history data files produced by ARPS.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!    12/12/2003: Rewrote based on arpscvt. Input file names will be
!    specified in arpsread.input.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nzsoil   ! Grid dimensions.
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER :: hinfmt,nhisfile_max,nhisfile,lengbf,nf,lenfil
  PARAMETER (nhisfile_max=200)
  CHARACTER (LEN=256) :: grdbasfn,hisfile(nhisfile_max)
!
!-----------------------------------------------------------------------
!  Include files:
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'timelvls.inc'
! INCLUDE 'indtflg.inc'
! INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!  Arrays to be read in:
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x(:)   ! The x-coord. of the physical and
                              ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y(:)   ! The y-coord. of the physical and
                              ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z(:)   ! The z-coord. of the computational
                              ! grid. Defined at w-point.

  REAL, ALLOCATABLE :: zp(:,:,:)  ! The height of the terrain.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)  ! The height of the terrain.

  REAL, ALLOCATABLE :: uprt   (:,:,:) ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:) ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:) ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:) ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:) ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:) ! Perturbation water vapor specific
                                      ! humidity (kg/kg)
!  REAL, ALLOCATABLE :: qc     (:,:,:) ! Cloud water mixing ratio (kg/kg)
!  REAL, ALLOCATABLE :: qr     (:,:,:) ! Rain water mixing ratio (kg/kg)
!  REAL, ALLOCATABLE :: qi     (:,:,:) ! Cloud ice mixing ratio (kg/kg)
!  REAL, ALLOCATABLE :: qs     (:,:,:) ! Snow mixing ratio (kg/kg)
!  REAL, ALLOCATABLE :: qh     (:,:,:) ! Hail mixing ratio (kg/kg)

  REAL, ALLOCATABLE :: qscalar (:,:,:,:) ! Scalar microphysics array

  REAL, ALLOCATABLE :: tke    (:,:,:) ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: kmh    (:,:,:) ! Horizontal turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv    (:,:,:) ! Vertical turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: ubar   (:,:,:) ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar   (:,:,:) ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar   (:,:,:) ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar  (:,:,:) ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar   (:,:,:) ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: rhobar (:,:,:) ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: qvbar  (:,:,:) ! Base state water vapor specific
                                      ! humidity (kg/kg)

  REAL, ALLOCATABLE :: u      (:,:,:) ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v      (:,:,:) ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w      (:,:,:) ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: qv     (:,:,:) ! Water vapor specific humidity (kg/kg)

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type
  REAL, ALLOCATABLE :: stypfrct(:,:,:)    ! Soil type
  INTEGER, ALLOCATABLE :: vegtyp(:,:)     ! Vegetation type
  REAL, ALLOCATABLE :: lai    (:,:)   ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns (:,:)   ! Surface roughness
  REAL, ALLOCATABLE :: veg    (:,:)   ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil (:,:,:,:) ! Soil temperature (K)
  REAL, ALLOCATABLE :: qsoil (:,:,:,:) ! Soil moisture (m**3/m**3)
  REAL, ALLOCATABLE :: wetcanp(:,:,:) ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)  ! Snow depth (m)

  REAL, ALLOCATABLE :: raing(:,:)     ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc(:,:)     ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:) ! precipitation rate (kg/(m**2*s))
                                      ! prcrate(1,1,1) = total precip. rate
                                      ! prcrate(1,1,2) = grid scale precip. rate
                                      ! prcrate(1,1,3) = cumulus precip. rate
                                      ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:)  ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw (:,:)    ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)    ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet (:,:) ! Net solar radiation, SWin - SWout
  REAL, ALLOCATABLE :: radlwin  (:,:) ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)    ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)    ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)    ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)    ! Surface moisture flux (kg/(m**2*s))

  REAL, ALLOCATABLE :: reflectivity(:,:,:)  ! simulated reflectivity (dBZ)
!
!-----------------------------------------------------------------------
! Temporary working arrays:
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
!
!-----------------------------------------------------------------------
! Misc. internal variables
!-----------------------------------------------------------------------
!
  INTEGER :: nchin

  INTEGER :: grdbas
  INTEGER :: i,j,k,nq,ireturn

  REAL :: time

  INTEGER :: nfile

  INTEGER :: length, ierr, istatus

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!  Get the names of the input data files.
!-----------------------------------------------------------------------
!
  CALL get_input_file_names(5,hinfmt,grdbasfn,hisfile,nhisfile)

  lenfil = len_trim(hisfile(1))

  CALL get_dims_from_data(hinfmt,hisfile(1)(1:lenfil),                    &
       nx,ny,nz,nzsoil,nstyps, ireturn)

  IF (nstyps <= 0) nstyps = 1
  nstyp = nstyps

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(4(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz ,'nzsoil=',nzsoil

  ALLOCATE(x      (nx))
  ALLOCATE(y      (ny))
  ALLOCATE(z      (nz))
  ALLOCATE(zp     (nx,ny,nz))
  ALLOCATE(zpsoil (nx,ny,nzsoil))

  ALLOCATE(uprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "uprt")

  ALLOCATE(vprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "vprt")

  ALLOCATE(wprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "wprt")

  ALLOCATE(ptprt  (nx,ny,nz))
  ALLOCATE(pprt   (nx,ny,nz))
  ALLOCATE(qvprt  (nx,ny,nz))
!  ALLOCATE(qc     (nx,ny,nz))
!  ALLOCATE(qr     (nx,ny,nz))
!  ALLOCATE(qi     (nx,ny,nz))
!  ALLOCATE(qs     (nx,ny,nz))
!  ALLOCATE(qh     (nx,ny,nz))

  ALLOCATE(qscalar(nx,ny,nz,nscalar))
  qscalar(:,:,:,:) = 0.0

  ALLOCATE(tke    (nx,ny,nz))
  ALLOCATE(kmh    (nx,ny,nz))
  ALLOCATE(kmv    (nx,ny,nz))
  ALLOCATE(ubar   (nx,ny,nz))
  ALLOCATE(vbar   (nx,ny,nz))
  ALLOCATE(wbar   (nx,ny,nz))
  ALLOCATE(ptbar  (nx,ny,nz))
  ALLOCATE(pbar   (nx,ny,nz))
  ALLOCATE(rhobar (nx,ny,nz))

  ALLOCATE(qvbar  (nx,ny,nz))
  ALLOCATE(u      (nx,ny,nz))
  ALLOCATE(v      (nx,ny,nz))
  ALLOCATE(w      (nx,ny,nz))
  ALLOCATE(qv     (nx,ny,nz))

  ALLOCATE(soiltyp (nx,ny,nstyps))
  ALLOCATE(stypfrct(nx,ny,nstyps))
  ALLOCATE(vegtyp (nx,ny))
  ALLOCATE(lai    (nx,ny))
  ALLOCATE(roufns (nx,ny))
  ALLOCATE(veg    (nx,ny))

  ALLOCATE(tsoil (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(qsoil (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(wetcanp(nx,ny,0:nstyps))
  ALLOCATE(snowdpth(nx,ny))

  ALLOCATE(raing(nx,ny))
  ALLOCATE(rainc(nx,ny))
  ALLOCATE(prcrate(nx,ny,4))

  ALLOCATE(radsw (nx,ny))
  ALLOCATE(rnflx (nx,ny))
  ALLOCATE(radswnet (nx,ny))
  ALLOCATE(radlwin (nx,ny))

  ALLOCATE(radfrc(nx,ny,nz))

  ALLOCATE(usflx (nx,ny))
  ALLOCATE(vsflx (nx,ny))
  ALLOCATE(ptsflx(nx,ny))
  ALLOCATE(qvsflx(nx,ny))
  ALLOCATE(tem1(nx,ny,nz))
  ALLOCATE(tem2(nx,ny,nz))
  ALLOCATE(tem3(nx,ny,nz))

  ALLOCATE(reflectivity(nx,ny,nz))

  x      =0.0
  y      =0.0
  z      =0.0
  zp     =0.0
  zpsoil =0.0

  uprt   =0.0
  vprt   =0.0
  wprt   =0.0
  ptprt  =0.0
  pprt   =0.0
  qvprt  =0.0
!  qc     =0.0
!  qr     =0.0
!  qi     =0.0
!  qs     =0.0
!  qh     =0.0
  tke    =0.0
  kmh    =0.0
  kmv    =0.0
  ubar   =0.0
  vbar   =0.0
  wbar   =0.0
  ptbar  =0.0
  pbar   =0.0
  rhobar =0.0
  qvbar  =0.0
  u      =0.0
  v      =0.0
  w      =0.0
  qv     =0.0

  soiltyp =0.0
  stypfrct=0.0
  vegtyp =0.0
  lai    =0.0
  roufns =0.0
  veg    =0.0

  tsoil  =0.0
  qsoil  =0.0
  wetcanp=0.0
  snowdpth=0.0

  raing=0.0
  rainc=0.0
  prcrate=0.0

  radfrc=0.0
  radsw =0.0
  rnflx =0.0
  radswnet = 0.0
  radlwin = 0.0

  usflx =0.0
  vsflx =0.0
  ptsflx=0.0
  qvsflx=0.0
  tem1=0.0
  tem2=0.0
  tem3=0.0
!
!-----------------------------------------------------------------------
!  Get the name of the input data set.
!-----------------------------------------------------------------------
!
  ldirnam=LEN(dirname)
  CALL strlnth( dirname , ldirnam)

  lengbf=len_trim(grdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', grdbasfn(1:lengbf)

  DO nfile = 1,nhisfile

    lenfil=len_trim(hisfile(nfile))
    WRITE(6,'(/a,a,a)')                                                 &
        ' Data set ', trim(hisfile(nfile)),' to be converted.'
!
!-----------------------------------------------------------------------
!  Read all input data arrays
!-----------------------------------------------------------------------

    CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                &
                 hinfmt, nchin,grdbasfn(1:lengbf),lengbf,               &
                 hisfile(nfile)(1:lenfil),lenfil,time,                  &
                 x,y,z,zp,zpsoil, uprt ,vprt ,wprt ,ptprt, pprt ,       &
                 qvprt, qscalar(:,:,:,:), tke,kmh,kmv,                  &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil, qsoil, wetcanp,snowdpth,                        &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 ireturn, tem1, tem2, tem3)

    curtim = time

    IF( ireturn /= 0 ) THEN
       WRITE(6,'(1x,a,i2,/1x,a)')                                       &
      'Data read was unsuccessful. ireturn =', ireturn,' Job stopped.'
       STOP
    END IF

!
!   Do whatever you want to do with the data here
!

  END DO

  STOP
END PROGRAM ARPSREAD
