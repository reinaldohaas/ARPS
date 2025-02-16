!This file contains all IO API for the netCDF format
!
! Author: Yunheng Wang (09/17/2004)
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE NETREAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netread(netid,packed,itime,grdbas,time,                      &
                 nx,ny,nz,nzsoil,nstyps, x, y, z, zp,zpsoil,            &
                 uprt, vprt, wprt, ptprt, pprt, qvprt,                  &
                 qscalar, tke,kmh,kmv,                                  &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1, ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read ARPS history data from NetCDF file.
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
  INCLUDE   'indtflg.inc'
  INCLUDE   'globcst.inc'
  INCLUDE   'grid.inc'          ! Grid & map parameters.
  INCLUDE   'phycst.inc'
  INCLUDE   'mp.inc'            ! mpi parameters.
!
!-----------------------------------------------------------------------
!
! Variable decalarations
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: netid
  INTEGER, INTENT(IN)  :: packed
  INTEGER, INTENT(IN)  :: itime
  INTEGER, INTENT(IN)  :: grdbas               ! Data read flag.
  INTEGER, INTENT(IN)  :: nx,ny,nz             ! Number of grid points in 3 directions
  INTEGER, INTENT(IN)  :: nzsoil               ! Number of grid points in the soil
  INTEGER, INTENT(IN)  :: nstyps               ! Number of soil type

  REAL,    INTENT(OUT) :: time                 ! Time in seconds of data read
                                               ! from "filename"
  REAL,    INTENT(OUT) :: x     (nx)           ! x-coord. of the physical and compu
                                               ! -tational grid. Defined at u-point(m).
  REAL,    INTENT(OUT) :: y     (ny)           ! y-coord. of the physical and compu
                                               ! -tational grid. Defined at v-point(m).
  REAL,    INTENT(OUT) :: z     (nz)           ! z-coord. of the computational grid.
                                               ! Defined at w-point on the staggered
                                               ! grid(m).
  REAL,    INTENT(OUT) :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                                               ! w-point of the staggered grid(m).
  REAL,    INTENT(OUT) :: zpsoil(nx,ny,nzsoil) ! Physical height coordinate defined at
                                               ! w-point of the soil (m)
  REAL,    INTENT(OUT) :: uprt  (nx,ny,nz)     ! Perturbation u-velocity (m/s)
  REAL,    INTENT(OUT) :: vprt  (nx,ny,nz)     ! Perturbation v-velocity (m/s)
  REAL,    INTENT(OUT) :: wprt  (nx,ny,nz)     ! Perturbation w-velocity (m/s)
  REAL,    INTENT(OUT) :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL,    INTENT(OUT) :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL,    INTENT(OUT) :: qvprt (nx,ny,nz)     ! Perturbation water vapor mixing
                                               ! ratio (kg/kg)

  REAL,    INTENT(OUT) :: qscalar    (nx,ny,nz,nscalar)

  REAL,    INTENT(OUT) :: tke  (nx,ny,nz)      ! Turbulent Kinetic Energy ((m/s)**2)
  REAL,    INTENT(OUT) :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                                               ! momentum. ( m**2/s )
  REAL,    INTENT(OUT) :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                                               ! momentum. ( m**2/s )

  REAL,    INTENT(INOUT) :: ubar  (nx,ny,nz)   ! Base state u-velocity (m/s)
  REAL,    INTENT(INOUT) :: vbar  (nx,ny,nz)   ! Base state v-velocity (m/s)
  REAL,    INTENT(INOUT) :: wbar  (nx,ny,nz)   ! Base state w-velocity (m/s)
  REAL,    INTENT(INOUT) :: ptbar (nx,ny,nz)   ! Base state potential temperature (K)
  REAL,    INTENT(INOUT) :: pbar  (nx,ny,nz)   ! Base state pressure (Pascal)
  REAL,    INTENT(INOUT) :: rhobar(nx,ny,nz)   ! Base state air density (kg/m**3)
  REAL,    INTENT(INOUT) :: qvbar (nx,ny,nz)   ! Base state water vapor mixing ratio

  INTEGER, INTENT(OUT) :: soiltyp (nx,ny,nstyps)         ! Soil type
  REAL,    INTENT(OUT) :: stypfrct(nx,ny,nstyps)         ! Soil type fraction
  INTEGER, INTENT(OUT) :: vegtyp (nx,ny)                 ! Vegetation type
  REAL,    INTENT(OUT) :: lai    (nx,ny)                 ! Leaf Area Index
  REAL,    INTENT(OUT) :: roufns (nx,ny)                 ! Surface roughness
  REAL,    INTENT(OUT) :: veg    (nx,ny)                 ! Vegetation fraction

  REAL,    INTENT(OUT) :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL,    INTENT(OUT) :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL,    INTENT(OUT) :: wetcanp(nx,ny,0:nstyps)        ! Canopy water amount
  REAL,    INTENT(OUT) :: snowdpth(nx,ny)                ! Snow depth (m)

  REAL,    INTENT(OUT) :: raing(nx,ny)         ! Grid supersaturation rain
  REAL,    INTENT(OUT) :: rainc(nx,ny)         ! Cumulus convective rain
  REAL,    INTENT(OUT) :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL,    INTENT(OUT) :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL,    INTENT(OUT) :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL,    INTENT(OUT) :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL,    INTENT(OUT) :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL,    INTENT(OUT) :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL,    INTENT(OUT) :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL,    INTENT(OUT) :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL,    INTENT(OUT) :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
  REAL,    INTENT(OUT) :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL,    INTENT(INOUT) :: tem1(nx,ny,nz)     ! Temporary work array

  INTEGER, INTENT(OUT) :: ireturn              ! Return status indicator
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,is,n
  INTEGER :: nxin,nyin,nzin,nzsoilin
  INTEGER :: bgrdin,bbasin,bvarin,bicein,btkein,btrbin
  INTEGER :: idummy,nstyps1
  INTEGER, SAVE :: nstypsin

  INTEGER :: nq
  INTEGER :: nqscalarin(nscalar)

  CHARACTER(LEN=20) :: varname

  REAL,    ALLOCATABLE :: invar2d (:,:)
  REAL,    ALLOCATABLE :: invar3dt(:,:,:)
  REAL,    ALLOCATABLE :: invar3du(:,:,:)
  REAL,    ALLOCATABLE :: invar3dv(:,:,:)
  REAL,    ALLOCATABLE :: invar3dw(:,:,:)
  REAL,    ALLOCATABLE :: invar4d (:,:,:,:)
  INTEGER, ALLOCATABLE :: invar2di(:,:)
  INTEGER, ALLOCATABLE :: invar3di(:,:,:)

  CHARACTER(LEN=4) :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  IF ( mp_opt /= 1 ) readsplit = 0  ! NO-MPI or is not initialized

!-----------------------------------------------------------------------
!
!  Read dimensions and global attributes
!
!-----------------------------------------------------------------------
!
  IF ( itime == 1 ) THEN

    CALL netreaddims(netid,nxin,nyin,nzin,nzsoilin,nstypsin,ireturn)
    !
    ! Data validation: dimensions
    !
    IF( nxin /= nx .OR. nyin /= ny .OR.                                 &
        nzin /= nz .OR. nzsoil /= nzsoil) THEN
      WRITE(6,'(1x,a)')                                                 &
         ' Dimensions in NETREAD inconsistent with data.'
      WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin, nzsoilin
      WRITE(6,'(1x,a,3I15)') ' Expected:  ', nx,   ny,   nz,   nzsoil
      WRITE(6,'(1x,a)')      ' Program aborted in NETREAD.'
      CALL arpsstop('arpstop called from binread nx-ny-nz read ',1)
    END IF

    IF (nstypsin > nstyps) THEN

      WRITE(6,'(/1x,3(a,I4),a/)') 'WARNING: nstyps in the data file is ',&
             nstypsin, ' which is larger than the decalared dimension ',&
             nstyps,' only ',nstyps, ' soil types will be extracted.'

    ELSE if (nstypsin < nstyps) THEN

      WRITE(6,'(/1x,a,I4,a,a,I4,a/)') 'WARNING: only ',nstypsin,        &
                     ' soil types are available inside the data file.', &
                     ' Because the decalared dimension is ',nstyps,     &
                     ' the extra soil types will be packed with zeros.'

    END IF

    IF (grdbas == 1) THEN
      CALL netreadatts(netid,grdbas,runname,nocmnt,cmnt,dx,dy,          &
                       year,month,day,hour,minute,second,thisdmp,tstop, &
                       mapproj,sclfct,trulat1,trulat2,trulon,latitud,   &
                       ctrlat,ctrlon,xgrdorg,ygrdorg,umove,vmove,       &
                       ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice,      &
                       rhosnow,rhogrpl,rhohail,alpharain,alphaice,      &
                       alphasnow,alphagrpl,alphahail,                   &
                       bgrdin,bbasin,bvarin,mstin,bicein,btrbin,        &
                       idummy,idummy,landin,totin,btkein,               &
                       prcin,radin,flxin,snowin,                        &
                       nscalarin,nqscalarin,ireturn)
    ELSE
      CALL netreadatts(netid,grdbas,runname,nocmnt,cmnt,dx,dy,          &
                       year,month,day,hour,minute,second,thisdmp,tstop, &
                       mapproj,sclfct,trulat1,trulat2,trulon,latitud,   &
                       ctrlat,ctrlon,xgrdorg,ygrdorg,umove,vmove,       &
                       ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice,      &
                       rhosnow,rhogrpl,rhohail,alpharain,alphaice,      &
                       alphasnow,alphagrpl,alphahail,                   &
                       grdin,basin,varin,mstin,icein,trbin,             &
                       sfcin,rainin,landin,totin,tkein,                 &
                       prcin,radin,flxin,snowin,                        &
                       nscalarin,nqscalarin,ireturn)

    END IF

  END IF
  snowcin = 0
  nstyps1 = MAX(1, MIN(nstypsin,nstyps))    ! nstyps:   Decalared dimension
                                            ! nstypsin: Dimension inside file
                                            ! nstyps1:  Dimension to be extracted
  ALLOCATE(invar2d(nx-1,ny-1),                       STAT = idummy)
  CALL check_alloc_status(idummy, "NETREAD:invar2d")
  ALLOCATE(invar3dt(nx-1,ny-1,nz-1),                 STAT = idummy)
  CALL check_alloc_status(idummy, "NETREAD:invar3dt")
  ALLOCATE(invar3du(nx,  ny-1,nz-1),                 STAT = idummy)
  CALL check_alloc_status(idummy, "NETREAD:invar3dt")
  ALLOCATE(invar3dv(nx-1,ny,  nz-1),                 STAT = idummy)
  CALL check_alloc_status(idummy, "NETREAD:invar3dt")
  ALLOCATE(invar3dw(nx-1,ny-1,MAX(nz,nzsoil,nstypsin+1)), STAT = idummy)
  CALL check_alloc_status(idummy, "NETREAD:invar3dt")
  ALLOCATE(invar4d (nx-1,ny-1,nzsoil,nstypsin+1),     STAT = idummy)
  CALL check_alloc_status(idummy, "NETREAD:invar4d")

  ALLOCATE(invar2di(nx-1,ny-1),                      STAT = idummy)
  CALL check_alloc_status(idummy, "NETREAD:invar2di")
  ALLOCATE(invar3di(nx-1,ny-1,nstypsin),             STAT = idummy)
  CALL check_alloc_status(idummy, "NETREAD:invar3di")
!
!-----------------------------------------------------------------------
!
! Readin current valid time
!
!-----------------------------------------------------------------------

  IF ( grdbas /= 1 ) THEN
    CALL netreadTime(netid,itime,'Time',time)
    WRITE(6,'(1x,/,a,f8.1,a,f8.3,a/)')   'To read data for time:',      &
                                   time,' secs = ',(time/60.),' mins.'
  ELSE

    WRITE(6,'(1x,/,a,/)') 'To read grid and base state data.'

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in x, y, z and zp arrays.
!
!----------------------------------------------------------------------
!
  IF( itime <= 1 .AND. (grdin == 1 .OR. grdbas == 1) ) THEN

    CALL netread1d(netid,packed,0,'x_stag',nx,x)
    CALL netread1d(netid,packed,0,'y_stag',ny,y)
    CALL netread1d(netid,packed,0,'z_stag',nz,z)

    CALL netread3d(netid,packed,0,'ZP',nx-1,ny-1,nz,invar3dw)
    DO k = 1, nz
      DO j = 1, ny-1
        DO i = 1, nx-1
          zp(i,j,k) = invar3dw(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(zp,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)

    CALL netread3d(netid,packed,0,'ZPSOIL',nx-1,ny-1,nzsoil,invar3dw)
    DO k = 1, nzsoil
      DO j = 1, ny-1
        DO i = 1, nx-1
          zpsoil(i,j,k) = invar3dw(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(zpsoil,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nzsoil,1,nzsoil)

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!----------------------------------------------------------------------
!
  IF(itime == 1 .AND. (basin == 1 .OR. grdbas == 1) ) THEN

    CALL netread3d(netid,packed,0,'UBAR',nx,ny-1,nz-1,invar3du)
    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx
          ubar(i,j,k) = invar3du(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(ubar,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)


    CALL netread3d(netid,packed,0,'VBAR',nx-1,ny,nz-1,invar3dv)
    DO k = 1, nz-1
      DO j = 1, ny
        DO i = 1, nx-1
          vbar(i,j,k) = invar3dv(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(vbar,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)

    CALL netread3d(netid,packed,0,'WBAR',nx-1,ny-1,nz,invar3dw)
    DO k = 1, nz
      DO j = 1, ny-1
        DO i = 1, nx-1
          wbar(i,j,k) = invar3dw(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(wbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)

    CALL netread3d(netid,packed,0,'PTBAR',nx-1,ny-1,nz-1,invar3dt)
    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          ptbar(i,j,k) = invar3dt(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(ptbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    CALL netread3d(netid,packed,0,'PBAR',nx-1,ny-1,nz-1,invar3dt)
    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          pbar(i,j,k) = invar3dt(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(pbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    IF (mstin == 1) THEN

      CALL netread3d(netid,packed,0,'QVBAR',nx-1,ny-1,nz-1,invar3dt)
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            qvbar(i,j,k) = invar3dt(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(qvbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    END IF

    IF (landin == 1) THEN

      CALL netread3di(netid,packed,0,'SOILTYP',nx-1,ny-1,nstypsin,invar3di)
      DO is = 1, nstyps1
        DO j = 1, ny-1
          DO i = 1, nx-1
            soiltyp(i,j,is) = invar3di(i,j,is)
          END DO
        END DO
      END DO
      CALL iedgfill(soiltyp,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nstyps,1,nstyps1)

      CALL netread3d(netid,packed,0,'STYPFRCT',nx-1,ny-1,nstypsin,invar3dw)
      DO is = 1, nstyps1
        DO j = 1, ny-1
          DO i = 1, nx-1
            stypfrct(i,j,is) = invar3dw(i,j,is)
          END DO
        END DO
      END DO
      CALL fix_stypfrct_nstyp(nx,ny,nstyps1,nstyps,stypfrct)

      CALL netread2di(netid,packed,0,'VEGTYP',nx-1,ny-1,invar2di)
      DO j = 1, ny-1
        DO i = 1, nx-1
          vegtyp(i,j) = invar2di(i,j)
        END DO
      END DO
      CALL iedgfill(vegtyp,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

      CALL netread2d(netid,packed,0,'LAI',nx-1,ny-1,invar2d)
      DO j = 1, ny-1
        DO i = 1, nx-1
          lai(i,j) = invar2d(i,j)
        END DO
      END DO
      CALL edgfill(lai,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

      CALL netread2d(netid,packed,0,'ROUFNS',nx-1,ny-1,invar2d)
      DO j = 1, ny-1
        DO i = 1, nx-1
          roufns(i,j) = invar2d(i,j)
        END DO
      END DO
      CALL edgfill(roufns,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

      CALL netread2d(netid,packed,0,'VEG',nx-1,ny-1,invar2d)
      DO j = 1, ny-1
        DO i = 1, nx-1
          veg(i,j) = invar2d(i,j)
        END DO
      END DO
      CALL edgfill(veg,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

    END IF

  END IF

  IF ( grdbas == 1 ) GOTO 4444

  IF (varin == 1) THEN

    IF (totin == 0) THEN

      CALL netread3d(netid,packed,itime,'UPRT',nx,ny-1,nz-1,invar3du)
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx
            uprt(i,j,k) = invar3du(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(uprt,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)

      CALL netread3d(netid,packed,itime,'VPRT',nx-1,ny,nz-1,invar3dv)
      DO k = 1, nz-1
        DO j = 1, ny
          DO i = 1, nx-1
            vprt(i,j,k) = invar3dv(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(vprt,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)

      CALL netread3d(netid,packed,itime,'WPRT',nx-1,ny-1,nz,invar3dw)
      DO k = 1, nz
        DO j = 1, ny-1
          DO i = 1, nx-1
            wprt(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(wprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)

      CALL netread3d(netid,packed,itime,'PTPRT',nx-1,ny-1,nz-1,invar3dt)
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            ptprt(i,j,k) = invar3dt(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(ptprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

      CALL netread3d(netid,packed,itime,'PPRT',nx-1,ny-1,nz-1,invar3dt)
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            pprt(i,j,k) = invar3dt(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(pprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    ELSE

      CALL netread3d(netid,packed,itime,'U',nx,ny-1,nz-1,invar3du)
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx
            uprt(i,j,k) = invar3du(i,j,k) - ubar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(uprt,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)

      CALL netread3d(netid,packed,itime,'V',nx-1,ny,nz-1,invar3dv)
      DO k = 1, nz-1
        DO j = 1, ny
          DO i = 1, nx-1
            vprt(i,j,k) = invar3dv(i,j,k) - vbar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(vprt,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)

      CALL netread3d(netid,packed,itime,'W',nx-1,ny-1,nz,invar3dw)
      DO k = 1, nz
        DO j = 1, ny-1
          DO i = 1, nx-1
            wprt(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(wprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)

      CALL netread3d(netid,packed,itime,'PT',nx-1,ny-1,nz-1,invar3dt)
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            ptprt(i,j,k) = invar3dt(i,j,k) - ptbar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(ptprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

      CALL netread3d(netid,packed,itime,'P',nx-1,ny-1,nz-1,invar3dt)
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            pprt(i,j,k) = invar3dt(i,j,k) - pbar(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(pprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    END IF

  END IF   ! varin

  IF (mstin == 1) THEN

    IF (totin == 0) THEN

      CALL netread3d(netid,packed,itime,'QVPRT',nx-1,ny-1,nz-1,invar3dt)
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            qvprt(i,j,k) = invar3dt(i,j,k)
          END DO
        END DO
      END DO

    ELSE

      CALL netread3d(netid,packed,itime,'QV',nx-1,ny-1,nz-1,invar3dt)
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            qvprt(i,j,k) = invar3dt(i,j,k) - qvbar(i,j,k)
          END DO
        END DO
      END DO

    END IF
    CALL edgfill(qvprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)


    DO nq = 1,nscalar
      IF (nqscalarin(nq) > 0)  THEN
        CALL netread3d(netid,packed,itime,upcase(qnames(nq)),nx-1,ny-1,nz-1,invar3dt)
        DO k = 1, nz-1
          DO j = 1, ny-1
            DO i = 1, nx-1
              qscalar(i,j,k,nq) = invar3dt(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      END IF
    END DO

    IF( rainin == 1 ) THEN

      CALL netread2d(netid,packed,itime,'RAING',nx-1,ny-1,invar2d)
      DO j = 1, ny-1
        DO i = 1, nx-1
          raing(i,j) = invar2d(i,j)
        END DO
      END DO
      CALL edgfill(raing,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

      CALL netread2d(netid,packed,itime,'RAINC',nx-1,ny-1,invar2d)
      DO j = 1, ny-1
        DO i = 1, nx-1
          rainc(i,j) = invar2d(i,j)
        END DO
      END DO
      CALL edgfill(rainc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

    END IF

    IF (prcin == 1) THEN

      DO n = 1,4
        WRITE(varname,'(a,I1)') 'PRCRATE',n

        CALL netread2d(netid,packed,itime,varname,nx-1,ny-1,invar2d)
        DO j = 1, ny-1
          DO i = 1, nx-1
            prcrate(i,j,n) = invar2d(i,j)
          END DO
        END DO

      END DO
      CALL edgfill(prcrate,1,nx,1,nx-1, 1,ny,1,ny-1, 1,4,1,4)

    END IF

  END IF

  IF( tkein == 1 ) THEN

    CALL netread3d(netid,packed,itime,'TKE',nx-1,ny-1,nz-1,invar3dt)
    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          tke(i,j,k) = invar3dt(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(tke,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

  END IF

  IF( trbin == 1 ) THEN

    CALL netread3d(netid,packed,itime,'KMH',nx-1,ny-1,nz-1,invar3dt)
    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          kmh(i,j,k) = invar3dt(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(kmh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    CALL netread3d(netid,packed,itime,'KMV',nx-1,ny-1,nz-1,invar3dt)
    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          kmv(i,j,k) = invar3dt(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(kmv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

  END IF

  IF (sfcin == 1) THEN
    !
    ! NOTE: Soil type dimensions
    !       nstyps:   Required in this run
    !       nstypsin: soil types inside the data file
    !       nstyps1:  = MIN(nstypsin, nstyps), soil types to be extracted
    !
    CALL netread4d(netid,packed,itime,'TSOIL',nx-1,ny-1,nzsoil,nstypsin+1,invar4d)
    DO is = 0,nstyps1
      DO k = 1, nzsoil
        DO j = 1, ny-1
          DO i = 1, nx-1
            tsoil(i,j,k,is) = invar4d(i,j,k,is+1)
          END DO
        END DO
      END DO
      CALL edgfill(tsoil(:,:,:,is),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nzsoil,1,nzsoil)
    END DO

    CALL netread4d(netid,packed,itime,'QSOIL',nx-1,ny-1,nzsoil,nstypsin+1,invar4d)
    DO is = 0,nstyps1
      DO k = 1, nzsoil
        DO j = 1, ny-1
          DO i = 1, nx-1
            qsoil(i,j,k,is) = invar4d(i,j,k,is+1)
          END DO
        END DO
      END DO
      CALL edgfill(tsoil(:,:,:,is),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nzsoil,1,nzsoil)
    END DO

    CALL netread3d(netid,packed,itime,'WETCANP',nx-1,ny-1,nstypsin+1,invar3dw)
    DO is = 0,nstyps1
      DO j = 1, ny-1
        DO i = 1, nx-1
          wetcanp(i,j,is) = invar3dw(i,j,is+1)
        END DO
      END DO
    END DO
    CALL edgfill(wetcanp,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nstyps1+1,1,nstyps1+1)

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstyps1,nstyps,tsoil,qsoil,wetcanp)

    IF (snowin == 1) THEN

      CALL netread2d(netid,packed,itime,'SNOWDPTH',nx-1,ny-1,invar2d)
      DO j = 1, ny-1
        DO i = 1, nx-1
          snowdpth(i,j) = invar2d(i,j)
        END DO
      END DO
      CALL edgfill(snowdpth,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

    END IF

  END IF

  IF (radin == 1) THEN

    CALL netread3d(netid,packed,itime,'RADFRC',nx-1,ny-1,nz-1,invar3dt)
    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          radfrc(i,j,k) = invar3dt(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(radfrc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    CALL netread2d(netid,packed,itime,'RADSW',nx-1,ny-1,invar2d)
    DO j = 1, ny-1
      DO i = 1, nx-1
        radsw(i,j) = invar2d(i,j)
      END DO
    END DO
    CALL edgfill(radsw,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

    CALL netread2d(netid,packed,itime,'RNFLX',nx-1,ny-1,invar2d)
    DO j = 1, ny-1
      DO i = 1, nx-1
        rnflx(i,j) = invar2d(i,j)
      END DO
    END DO
    CALL edgfill(rnflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

    CALL netread2d(netid,packed,itime,'RADSWNET',nx-1,ny-1,invar2d)
    DO j = 1, ny-1
      DO i = 1, nx-1
        radswnet(i,j) = invar2d(i,j)
      END DO
    END DO
    CALL edgfill(radswnet,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

    CALL netread2d(netid,packed,itime,'RADLWIN',nx-1,ny-1,invar2d)
    DO j = 1, ny-1
      DO i = 1, nx-1
        radlwin(i,j) = invar2d(i,j)
      END DO
    END DO
    CALL edgfill(radlwin,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

  END IF

  IF (flxin == 1) THEN

    CALL netread2d(netid,packed,itime,'USFLX',nx-1,ny-1,invar2d)
    DO j = 1, ny-1
      DO i = 1, nx-1
        usflx(i,j) = invar2d(i,j)
      END DO
    END DO
    CALL edgfill(usflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

    CALL netread2d(netid,packed,itime,'VSFLX',nx-1,ny-1,invar2d)
    DO j = 1, ny-1
      DO i = 1, nx-1
        vsflx(i,j) = invar2d(i,j)
      END DO
    END DO
    CALL edgfill(vsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

    CALL netread2d(netid,packed,itime,'PTSFLX',nx-1,ny-1,invar2d)
    DO j = 1, ny-1
      DO i = 1, nx-1
        ptsflx(i,j) = invar2d(i,j)
      END DO
    END DO
    CALL edgfill(ptsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

    CALL netread2d(netid,packed,itime,'QVSFLX',nx-1,ny-1,invar2d)
    DO j = 1, ny-1
      DO i = 1, nx-1
        qvsflx(i,j) = invar2d(i,j)
      END DO
    END DO
    CALL edgfill(qvsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)

  END IF

!-----------------------------------------------------------------------
!
! Clear memory and return
!
!-----------------------------------------------------------------------

  4444 CONTINUE

  DEALLOCATE(invar2d, invar2di, invar4d)
  DEALLOCATE(invar3dt,invar3du, invar3dv, invar3dw, invar3di)

  RETURN
END SUBROUTINE netread
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE NETREADSPLIT               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netreadsplit(netid,packed,itime,grdbas,time,                 &
                 nx,ny,nz,nzsoil,nstyps, x, y, z, zp,zpsoil,            &
                 uprt, vprt, wprt, ptprt, pprt, qvprt,                  &
                 qscalar, tke, kmh, kmv,                                &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1, ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read ARPS history data from NetCDF file.
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
  INCLUDE   'indtflg.inc'
  INCLUDE   'globcst.inc'
  INCLUDE   'phycst.inc'
  INCLUDE   'grid.inc'          ! Grid & map parameters.
  INCLUDE   'mp.inc'            ! mpi parameters.
!
!-----------------------------------------------------------------------
!
! Variable decalarations
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: netid
  INTEGER, INTENT(IN)  :: packed
  INTEGER, INTENT(IN)  :: itime
  INTEGER, INTENT(IN)  :: grdbas               ! Data read flag.
  INTEGER, INTENT(IN)  :: nx,ny,nz             ! Number of grid points in 3 directions
  INTEGER, INTENT(IN)  :: nzsoil               ! Number of grid points in the soil
  INTEGER, INTENT(IN)  :: nstyps               ! Number of soil type

  REAL,    INTENT(OUT) :: time                 ! Time in seconds of data read
                                               ! from "filename"
  REAL,    INTENT(OUT) :: x     (nx)           ! x-coord. of the physical and compu
                                               ! -tational grid. Defined at u-point(m).
  REAL,    INTENT(OUT) :: y     (ny)           ! y-coord. of the physical and compu
                                               ! -tational grid. Defined at v-point(m).
  REAL,    INTENT(OUT) :: z     (nz)           ! z-coord. of the computational grid.
                                               ! Defined at w-point on the staggered
                                               ! grid(m).
  REAL,    INTENT(OUT) :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                                               ! w-point of the staggered grid(m).
  REAL,    INTENT(OUT) :: zpsoil(nx,ny,nzsoil) ! Physical height coordinate defined at
                                               ! w-point of the soil (m)
  REAL,    INTENT(OUT) :: uprt  (nx,ny,nz)     ! Perturbation u-velocity (m/s)
  REAL,    INTENT(OUT) :: vprt  (nx,ny,nz)     ! Perturbation v-velocity (m/s)
  REAL,    INTENT(OUT) :: wprt  (nx,ny,nz)     ! Perturbation w-velocity (m/s)
  REAL,    INTENT(OUT) :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL,    INTENT(OUT) :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL,    INTENT(OUT) :: qvprt (nx,ny,nz)     ! Perturbation water vapor mixing
                                               ! ratio (kg/kg)

  REAL,    INTENT(OUT) :: qscalar    (nx,ny,nz,nscalar)

  REAL,    INTENT(OUT) :: tke  (nx,ny,nz)      ! Turbulent Kinetic Energy ((m/s)**2)
  REAL,    INTENT(OUT) :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                                               ! momentum. ( m**2/s )
  REAL,    INTENT(OUT) :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                                               ! momentum. ( m**2/s )

  REAL,    INTENT(INOUT) :: ubar  (nx,ny,nz)   ! Base state u-velocity (m/s)
  REAL,    INTENT(INOUT) :: vbar  (nx,ny,nz)   ! Base state v-velocity (m/s)
  REAL,    INTENT(INOUT) :: wbar  (nx,ny,nz)   ! Base state w-velocity (m/s)
  REAL,    INTENT(INOUT) :: ptbar (nx,ny,nz)   ! Base state potential temperature (K)
  REAL,    INTENT(INOUT) :: pbar  (nx,ny,nz)   ! Base state pressure (Pascal)
  REAL,    INTENT(INOUT) :: rhobar(nx,ny,nz)   ! Base state air density (kg/m**3)
  REAL,    INTENT(INOUT) :: qvbar (nx,ny,nz)   ! Base state water vapor mixing ratio

  INTEGER, INTENT(OUT) :: soiltyp (nx,ny,nstyps)         ! Soil type
  REAL,    INTENT(OUT) :: stypfrct(nx,ny,nstyps)         ! Soil type fraction
  INTEGER, INTENT(OUT) :: vegtyp (nx,ny)                 ! Vegetation type
  REAL,    INTENT(OUT) :: lai    (nx,ny)                 ! Leaf Area Index
  REAL,    INTENT(OUT) :: roufns (nx,ny)                 ! Surface roughness
  REAL,    INTENT(OUT) :: veg    (nx,ny)                 ! Vegetation fraction

  REAL,    INTENT(OUT) :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL,    INTENT(OUT) :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL,    INTENT(OUT) :: wetcanp(nx,ny,0:nstyps)        ! Canopy water amount
  REAL,    INTENT(OUT) :: snowdpth(nx,ny)                ! Snow depth (m)

  REAL,    INTENT(OUT) :: raing(nx,ny)         ! Grid supersaturation rain
  REAL,    INTENT(OUT) :: rainc(nx,ny)         ! Cumulus convective rain
  REAL,    INTENT(OUT) :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL,    INTENT(OUT) :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL,    INTENT(OUT) :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL,    INTENT(OUT) :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL,    INTENT(OUT) :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL,    INTENT(OUT) :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL,    INTENT(OUT) :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL,    INTENT(OUT) :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL,    INTENT(OUT) :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
  REAL,    INTENT(OUT) :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL,  INTENT(INOUT) :: tem1(nx,ny,nz)     ! Temporary work array

  INTEGER, INTENT(OUT) :: ireturn              ! Return status indicator
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,is,n
  INTEGER :: nxin,nyin,nzin,nzsoilin
  INTEGER :: bgrdin,bbasin,bvarin,bicein,btkein,btrbin
  INTEGER :: idummy,nstyps1
  INTEGER :: nxlg, nylg
  INTEGER, SAVE :: nstypsin

  INTEGER :: nq
  INTEGER :: nqscalarin(nscalar)

  CHARACTER(LEN=20) :: varname

  REAL,    ALLOCATABLE :: invar1d (:)       ! Used to extract data from NetCDF file
  REAL,    ALLOCATABLE :: invar3du(:,:,:)   ! U, UBAR
  REAL,    ALLOCATABLE :: invar3dv(:,:,:)   ! V, VBAR
  REAL,    ALLOCATABLE :: invar3dw(:,:,:)   ! W, WBAR, ZP,
                                            ! ZPSOIL, stypfrct, wetcanp
                                            ! and other 2D/3D nostag data
  REAL,    ALLOCATABLE :: invar4d (:,:,:,:) ! tsoil, qsoil
  INTEGER, ALLOCATABLE :: invar3di(:,:,:)   ! soiltyp, vegtyp

  REAL,    ALLOCATABLE :: var3d (:,:,:)     ! used to split those data
  INTEGER, ALLOCATABLE :: var3di(:,:,:)
  REAL,    ALLOCATABLE :: var4d (:,:,:,:)

  CHARACTER(LEN=4)  :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

!-----------------------------------------------------------------------
!
!  Read dimensions and global attributes
!
!-----------------------------------------------------------------------
!
  IF ( itime == 1 .AND. myproc == 0  ) THEN

    CALL netreaddims(netid,nxin,nyin,nzin,nzsoilin,nstypsin,ireturn)
    !
    ! Data validation: dimensions
    !
    IF( nxin /= nxlg .OR. nyin   /= nylg .OR.                           &
        nzin /= nz   .OR. nzsoil /= nzsoil) THEN
      WRITE(6,'(1x,a)')                                                 &
         ' Dimensions in NETREAD inconsistent with data.'
      WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin, nzsoilin
      WRITE(6,'(1x,a,3I15)') ' Expected:  ', nxlg, nylg,  nz,   nzsoil
      WRITE(6,'(1x,a)')      ' Program aborted in NETREAD.'
      CALL arpsstop('arpstop called from binread nx-ny-nz read ',1)
    END IF

    IF (nstypsin > nstyps) THEN

      WRITE(6,'(/1x,3(a,I4),a/)') 'WARNING: nstyps in the data file is ',&
             nstypsin, ' which is larger than the decalared dimension ',&
             nstyps,' only ',nstyps, ' soil types will be extracted.'

    ELSE if (nstypsin < nstyps) THEN

      WRITE(6,'(/1x,a,I4,a,a,I4,a/)') 'WARNING: only ',nstypsin,        &
                     ' soil types are available inside the data file.', &
                     ' Because the decalared dimension is ',nstyps,     &
                     ' the extra soil types will be packed with zeros.'

    END IF

    IF (grdbas == 1) THEN
      CALL netreadatts(netid,grdbas,runname,nocmnt,cmnt,dx,dy,          &
                       year,month,day,hour,minute,second,thisdmp,tstop, &
                       mapproj,sclfct,trulat1,trulat2,trulon,latitud,   &
                       ctrlat,ctrlon,xgrdorg,ygrdorg,umove,vmove,       &
                       ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice,      &
                       rhosnow,rhogrpl,rhohail,alpharain,alphaice,      &
                       alphasnow,alphagrpl,alphahail,                   &
                       bgrdin,bbasin,bvarin,mstin,bicein,btrbin,        &
                       idummy,idummy,landin,totin,btkein,               &
                       prcin,radin,flxin,snowin,                        &
                       nscalarin,nqscalarin,ireturn)
    ELSE
      CALL netreadatts(netid,grdbas,runname,nocmnt,cmnt,dx,dy,          &
                       year,month,day,hour,minute,second,thisdmp,tstop, &
                       mapproj,sclfct,trulat1,trulat2,trulon,latitud,   &
                       ctrlat,ctrlon,xgrdorg,ygrdorg,umove,vmove,       &
                       ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice,      &
                       rhosnow,rhogrpl,rhohail,alpharain,alphaice,      &
                       alphasnow,alphagrpl,alphahail,                   &
                       grdin,basin,varin,mstin,icein,trbin,             &
                       sfcin,rainin,landin,totin,tkein,                 &
                       prcin,radin,flxin,snowin,                        &
                       nscalarin,nqscalarin,ireturn)

    END IF

  END IF

  IF (itime == 1) THEN
    CALL mpupdatec(runname, 40)
    CALL mpupdatei(mstin,1)
    CALL mpupdatei(landin,1)
    CALL mpupdatei(totin,1)
    CALL mpupdatei(mapproj,1)
    CALL mpupdatei(month,1)
    CALL mpupdatei(day,1)
    CALL mpupdatei(year,1)
    CALL mpupdatei(hour,1)
    CALL mpupdatei(minute,1)
    CALL mpupdatei(second,1)
    IF(grdbas == 1) THEN
      CALL mpupdatei(bgrdin,1)
      CALL mpupdatei(bbasin,1)
      CALL mpupdatei(bvarin,1)
      CALL mpupdatei(btrbin,1)
      CALL mpupdatei(btkein,1)
    ELSE
      CALL mpupdatei(grdin,1)
      CALL mpupdatei(basin,1)
      CALL mpupdatei(varin,1)
      CALL mpupdatei(trbin,1)
      CALL mpupdatei(tkein,1)
      CALL mpupdatei(icein,1)
      CALL mpupdatei(sfcin,1)
      CALL mpupdatei(rainin,1)

      CALL mpupdatei(nscalarin,1)
      CALL mpupdatei(nqscalarin,nscalar)
    END IF

    CALL mpupdater(umove,1)
    CALL mpupdater(vmove,1)
    CALL mpupdater(xgrdorg,1)
    CALL mpupdater(ygrdorg,1)
    CALL mpupdater(trulat1,1)
    CALL mpupdater(trulat2,1)
    CALL mpupdater(trulon,1)
    CALL mpupdater(sclfct,1)
    CALL mpupdater(tstop,1)
    CALL mpupdater(thisdmp,1)
    CALL mpupdater(latitud,1)
    CALL mpupdater(ctrlat,1)
    CALL mpupdater(ctrlon,1)

    CALL mpupdater(ntcloud,1)
    CALL mpupdater(n0rain,1)
    CALL mpupdater(n0snow,1)
    CALL mpupdater(n0grpl,1)
    CALL mpupdater(n0hail,1)
    CALL mpupdater(rhoice,1)
    CALL mpupdater(rhosnow,1)
    CALL mpupdater(rhogrpl,1)
    CALL mpupdater(rhohail,1)
    CALL mpupdater(alpharain,1)
    CALL mpupdater(alphaice,1)
    CALL mpupdater(alphasnow,1)
    CALL mpupdater(alphagrpl,1)
    CALL mpupdater(alphahail,1)

    IF(totin /= 0) THEN
      CALL mpupdatei(prcin,1)
      CALL mpupdatei(radin,1)
      CALL mpupdatei(flxin,1)
      CALL mpupdatei(snowin,1)
    END IF

    CALL mpupdatei(nstypsin,1)
  END IF

  nstyps1 = MAX(1, MIN(nstypsin,nstyps))    ! nstyps:   Decalared dimension
                                            ! nstypsin: Dimension inside file
                                            ! nstyps1:  Dimension to be extracted
  snowcin = 0

  ALLOCATE(invar1d (MAX(nxlg,nylg)),                         STAT = idummy)
  CALL check_alloc_status(idummy, "NETREADSPLIT:invar1d")
  ALLOCATE(invar3du(nxlg,nylg-1,nz-1),                       STAT = idummy)
  CALL check_alloc_status(idummy, "NETREADSPLIT:invar3du")
  ALLOCATE(invar3dv(nxlg-1,nylg,nz-1),                       STAT = idummy)
  CALL check_alloc_status(idummy, "NETREADSPLIT:invar3dv")
  ALLOCATE(invar3dw(nxlg-1,nylg-1,MAX(nz,nzsoil,nstypsin+1)), STAT = idummy)
  CALL check_alloc_status(idummy, "NETREADSPLIT:invar3dw")
  ALLOCATE(invar4d (nxlg-1,nylg-1,nzsoil,nstypsin+1),        STAT = idummy)
  CALL check_alloc_status(idummy, "NETREADSPLIT:invar4d")
  ALLOCATE(invar3di(nxlg-1,nylg-1,nstypsin),                 STAT = idummy)
  CALL check_alloc_status(idummy, "NETREADSPLIT:invar3di")

  ALLOCATE(var3d (nxlg,nylg,MAX(nz,nzsoil,nstyps1+1)),   STAT = idummy)
  CALL check_alloc_status(idummy, "NETREADSPLIT:var3d")
  ALLOCATE(var3di(nxlg,nylg,nstyps1),                    STAT = idummy)
  CALL check_alloc_status(idummy, "NETREADSPLIT:var3di")
  ALLOCATE(var4d(nxlg,nylg,nzsoil,nstyps1+1),            STAT = idummy)
  CALL check_alloc_status(idummy, "NETREADSPLIT:var4d")
!
!-----------------------------------------------------------------------
!
! Readin current valid time
!
!-----------------------------------------------------------------------

  IF ( grdbas /= 1 ) THEN

    IF (myproc == 0 ) THEN
      CALL netreadTime(netid,itime,'Time',time)
      WRITE(6,'(1x,/,a,f8.1,a,f8.3,a/)')   'To read data for time:',    &
                                     time,' secs = ',(time/60.),' mins.'
    END IF
    CALL mpupdater(time,1)

  ELSE

    IF (myproc == 0)     &
      WRITE(6,'(1x,/,a,/)') 'To read grid and base state data.'

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in x, y, z and zp arrays.
!
!----------------------------------------------------------------------
!
  IF( itime == 1 .AND. (grdin == 1 .OR. grdbas == 1) ) THEN

    IF ( myproc == 0 ) THEN
      CALL netread1d(netid,packed,0,'x_stag',nxlg,invar1d)
    END IF
    CALL mpisplit1dx(invar1d,nx,x)

    IF ( myproc == 0 ) THEN
      CALL netread1d(netid,packed,0,'y_stag',nylg,invar1d)
    END IF
    CALL mpisplit1dy(invar1d,ny,y)

    IF ( myproc == 0 ) THEN
      CALL netread1d(netid,packed,0,'z_stag',nz,z)
    END IF
    CALL mpupdater(z,nz)

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,0,'ZP',nxlg-1,nylg-1,nz,invar3dw)
      DO k = 1,nz
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,zp)

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,0,'ZPSOIL',nxlg-1,nylg-1,nzsoil,invar3dw)
      DO k = 1,nzsoil
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nzsoil,1,nzsoil)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nzsoil,zpsoil)

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!----------------------------------------------------------------------
!
  IF(itime == 1 .AND. (basin == 1 .OR. grdbas == 1) ) THEN

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,0,'UBAR',nxlg,nylg-1,nz-1,invar3du)
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg
            var3d(i,j,k) = invar3du(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg,1,nylg,1,nylg-1,1,nz,1,nz-1)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,ubar)

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,0,'VBAR',nxlg-1,nylg,nz-1,invar3dv)
      DO k = 1,nz-1
        DO j = 1,nylg
          DO i = 1,nxlg-1
            var3d(i,j,k) = invar3dv(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg,1,nz,1,nz-1)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,vbar)

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,0,'WBAR',nxlg-1,nylg-1,nz,invar3dw)
      DO k = 1,nz
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,wbar)

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,0,'PTBAR',nxlg-1,nylg-1,nz-1,invar3dw)
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,ptbar)

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,0,'PBAR',nxlg-1,nylg-1,nz-1,invar3dw)
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,pbar)

    IF (mstin == 1) THEN

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,0,'QVBAR',nxlg-1,nylg-1,nz-1,invar3dw)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF

      CALL mpisplit3d(var3d,nx,ny,nz,qvbar)

    END IF

    IF (landin == 1) THEN

      IF ( myproc == 0  ) THEN
        CALL netread3di(netid,packed,0,'SOILTYP',nxlg-1,nylg-1,nstypsin,invar3di)
        DO is = 1,nstyps1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3di(i,j,is) = invar3di(i,j,is)
            END DO
          END DO
        END DO
        CALL iedgfill(var3di,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nstyps1,1,nstyps1)
      END IF
      CALL mpisplit3di(invar3di,nx,ny,nstyps1,soiltyp)

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,0,'STYPFRCT',nxlg-1,nylg-1,nstypsin,invar3dw)
        DO is = 1,nstyps1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,is) = invar3dw(i,j,is)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nstyps1,1,nstyps1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nstyps1,stypfrct)

      CALL fix_stypfrct_nstyp(nx,ny,nstyps1,nstyps,stypfrct)

      IF ( myproc == 0  ) THEN
        CALL netread2di(netid,packed,0,'VEGTYP',nxlg-1,nylg-1,invar3di)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3di(i,j,1) = invar3di(i,j,1)
          END DO
        END DO
        CALL iedgfill(var3di,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2di(var3di,nx,ny,vegtyp)

      IF ( myproc == 0  ) THEN
        CALL netread2d(netid,packed,0,'LAI',nxlg-1,nylg-1,invar3dw)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = invar3dw(i,j,1)
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var3d,nx,ny,lai)

      IF ( myproc == 0  ) THEN
        CALL netread2d(netid,packed,0,'ROUFNS',nxlg-1,nylg-1,invar3dw)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = invar3dw(i,j,1)
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var3d,nx,ny,roufns)

      IF ( myproc == 0  ) THEN
        CALL netread2d(netid,packed,0,'VEG',nxlg-1,nylg-1,invar3dw)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = invar3dw(i,j,1)
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var3d,nx,ny,veg)

    END IF

  END IF

  IF ( grdbas == 1 ) GOTO 4444

  IF (varin == 1) THEN

    IF (totin == 0) THEN

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'UPRT',nxlg,nylg-1,nz-1,invar3du)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg
              var3d(i,j,k) = invar3du(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,uprt)

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'VPRT',nxlg-1,nylg,nz-1,invar3dv)
        DO k = 1,nz-1
          DO j = 1,nylg
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dv(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,vprt)

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'WPRT',nxlg-1,nylg-1,nz,invar3dw)
        DO k = 1,nz
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,wprt)

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'PTPRT',nxlg-1,nylg-1,nz-1,invar3dw)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,ptprt)

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'PPRT',nxlg-1,nylg-1,nz-1,invar3dw)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,pprt)

    ELSE

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'U',nxlg,nylg-1,nz-1,invar3du)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg
              var3d(i,j,k) = invar3du(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,uprt)
      DO k = 1,nz
        DO j = 1,ny
          DO i = 1,nx
            uprt(i,j,k) = uprt(i,j,k) - ubar(i,j,k)
          END DO
        END DO
      END DO

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'V',nxlg-1,nylg,nz-1,invar3dv)
        DO k = 1,nz-1
          DO j = 1,nylg
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dv(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,vprt)
      DO k = 1,nz
        DO j = 1,ny
          DO i = 1,nx
            vprt(i,j,k) = vprt(i,j,k) - vbar(i,j,k)
          END DO
        END DO
      END DO

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'W',nxlg-1,nylg-1,nz,invar3dw)
        DO k = 1,nz
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,wprt)

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'PT',nxlg-1,nylg-1,nz-1,invar3dw)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,ptprt)
      ptprt(:,:,:) = ptprt(:,:,:) - ptbar(:,:,:)

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'P',nxlg-1,nylg-1,nz-1,invar3dw)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,pprt)
      pprt(:,:,:) = pprt(:,:,:) - pbar(:,:,:)
    END IF

  END IF   ! varin

  IF (mstin == 1) THEN

    IF (totin == 0) THEN

      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'QVPRT',nxlg-1,nylg-1,nz-1,invar3dw)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,qvprt)
    ELSE
      IF ( myproc == 0  ) THEN
        CALL netread3d(netid,packed,itime,'QV',nxlg-1,nylg-1,nz-1,invar3dw)
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = invar3dw(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
      END IF
      CALL mpisplit3d(var3d,nx,ny,nz,qvprt)
      qvprt(:,:,:) = qvprt(:,:,:) - qvbar(:,:,:)
    END IF

    DO nq = 1,nscalar
      IF (nqscalarin(nq) > 0 ) THEN
        IF ( myproc == 0  ) THEN
          CALL netread3d(netid,packed,itime,upcase(qnames(nq)),nxlg-1,nylg-1,nz-1,invar3dw)
          DO k = 1,nz-1
            DO j = 1,nylg-1
              DO i = 1,nxlg-1
                var3d(i,j,k) = invar3dw(i,j,k)
              END DO
            END DO
          END DO
          CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
        END IF
        CALL mpisplit3d(var3d,nx,ny,nz,qscalar(:,:,:,nq))
      END IF
    END DO

    IF( rainin == 1 ) THEN

      IF ( myproc == 0  ) THEN
        CALL netread2d(netid,packed,itime,'RAING',nxlg-1,nylg-1,invar3dw)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = invar3dw(i,j,1)
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var3d,nx,ny,raing)

      IF ( myproc == 0  ) THEN
        CALL netread2d(netid,packed,itime,'RAINC',nxlg-1,nylg-1,invar3dw)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = invar3dw(i,j,1)
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var3d,nx,ny,rainc)

    END IF

    IF (prcin == 1) THEN

      DO n = 1,4
        WRITE(varname,'(a,I1)') 'PRCRATE',n

        IF ( myproc == 0  ) THEN
          CALL netread2d(netid,packed,itime,varname,nxlg-1,nylg-1,invar3dw)
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,1) = invar3dw(i,j,1)
            END DO
          END DO
          CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
        END IF
        CALL mpisplit2d(var3d,nx,ny,prcrate(:,:,n))
      END DO

    END IF

  END IF

  IF( tkein == 1 ) THEN

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,itime,'TKE',nxlg-1,nylg-1,nz-1,invar3dw)
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,tke)

  END IF

  IF( trbin == 1 ) THEN

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,itime,'KMH',nxlg-1,nylg-1,nz-1,invar3dw)
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,kmh)

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,itime,'KMV',nxlg-1,nylg-1,nz-1,invar3dw)
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,kmv)

  END IF

  IF (sfcin == 1) THEN
    !
    ! NOTE: Soil type dimensions
    !       nstyps:   Required in this run
    !       nstypsin: soil types inside the data file
    !       nstyps1:  = MIN(nstypsin, nstyps), soil types to be extracted
    !
    IF ( myproc == 0  ) THEN
      CALL netread4d(netid,packed,itime,'TSOIL',nxlg-1,nylg-1,nzsoil,nstypsin+1,invar4d)
      DO is = 1, nstyps1+1
        DO k = 1,nzsoil
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var4d(i,j,k,is) = invar4d(i,j,k,is)
            END DO
          END DO
        END DO
        CALL edgfill(var4d(:,:,:,is),1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nzsoil,1,nzsoil)
      END DO
    END IF
    CALL mpisplit4d(var4d,nx,ny,nzsoil,nstyps1+1,tsoil)

    IF ( myproc == 0  ) THEN
      CALL netread4d(netid,packed,itime,'QSOIL',nxlg-1,nylg-1,nzsoil,nstypsin+1,invar4d)
      DO is = 1, nstyps1+1
        DO k = 1,nzsoil
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var4d(i,j,k,is) = invar4d(i,j,k,is)
            END DO
          END DO
        END DO
        CALL edgfill(var4d(:,:,:,is),1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nzsoil,1,nzsoil)
      END DO
    END IF
    CALL mpisplit4d(var4d,nx,ny,nzsoil,nstyps1+1,qsoil)

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,itime,'WETCANP',nxlg-1,nylg-1,nstypsin+1,invar3dw)
      DO is = 1, nstyps1+1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,is) = invar3dw(i,j,is)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nstyps1+1,1,nstyps1+1)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nstyps1+1,wetcanp)

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstyps1,nstyps,tsoil,qsoil,wetcanp)

    IF (snowin == 1) THEN

      IF ( myproc == 0  ) THEN
        CALL netread2d(netid,packed,itime,'SNOWDPTH',nxlg-1,nylg-1,invar3dw)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = invar3dw(i,j,1)
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var3d,nx,ny,snowdpth)

    END IF

  END IF

  IF (radin == 1) THEN

    IF ( myproc == 0  ) THEN
      CALL netread3d(netid,packed,itime,'RADFRC',nxlg-1,nylg-1,nz-1,invar3dw)
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,k) = invar3dw(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nz,1,nz-1)
    END IF
    CALL mpisplit3d(var3d,nx,ny,nz,radfrc)

    IF ( myproc == 0  ) THEN
      CALL netread2d(netid,packed,itime,'RADSW',nxlg-1,nylg-1,invar3dw)
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d(i,j,1) = invar3dw(i,j,1)
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
    END IF
    CALL mpisplit2d(var3d,nx,ny,radsw)

    IF ( myproc == 0  ) THEN
      CALL netread2d(netid,packed,itime,'RNFLX',nxlg-1,nylg-1,invar3dw)
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d(i,j,1) = invar3dw(i,j,1)
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
    END IF
    CALL mpisplit2d(var3d,nx,ny,rnflx)

    IF ( myproc == 0  ) THEN
      CALL netread2d(netid,packed,itime,'RADSWNET',nxlg-1,nylg-1,invar3dw)
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d(i,j,1) = invar3dw(i,j,1)
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
    END IF
    CALL mpisplit2d(var3d,nx,ny,radswnet)

    IF ( myproc == 0  ) THEN
      CALL netread2d(netid,packed,itime,'RADLWIN',nxlg-1,nylg-1,invar3dw)
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d(i,j,1) = invar3dw(i,j,1)
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
    END IF
    CALL mpisplit2d(var3d,nx,ny,radlwin)

  END IF

  IF (flxin == 1) THEN

    IF ( myproc == 0  ) THEN
      CALL netread2d(netid,packed,itime,'USFLX',nxlg-1,nylg-1,invar3dw)
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d(i,j,1) = invar3dw(i,j,1)
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
    END IF
    CALL mpisplit2d(var3d,nx,ny,usflx)

    IF ( myproc == 0  ) THEN
      CALL netread2d(netid,packed,itime,'VSFLX',nxlg-1,nylg-1,invar3dw)
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d(i,j,1) = invar3dw(i,j,1)
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
    END IF
    CALL mpisplit2d(var3d,nx,ny,vsflx)

    IF ( myproc == 0  ) THEN
      CALL netread2d(netid,packed,itime,'PTSFLX',nxlg-1,nylg-1,invar3dw)
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d(i,j,1) = invar3dw(i,j,1)
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
    END IF
    CALL mpisplit2d(var3d,nx,ny,ptsflx)

    IF ( myproc == 0  ) THEN
      CALL netread2d(netid,packed,itime,'QVSFLX',nxlg-1,nylg-1,invar3dw)
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d(i,j,1) = invar3dw(i,j,1)
        END DO
      END DO
      CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
    END IF
    CALL mpisplit2d(var3d,nx,ny,qvsflx)

  END IF

!-----------------------------------------------------------------------
!
! Clear memory and return
!
!-----------------------------------------------------------------------

  4444 CONTINUE

  DEALLOCATE(invar1d)
  DEALLOCATE(invar3du, invar3dv, invar3dw, invar4d)
  DEALLOCATE(var3di,invar3di)
  DEALLOCATE(var3d,  var4d)

  RETURN
END SUBROUTINE netreadsplit
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE NETDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE netdump(netid,itime,packed,nx,ny,nz,nzsoil,nstyps,grdbas,    &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 kmh,kmv,ubar,vbar,ptbar,pbar,rhobar,qvbar,             &
                 x,y,z,zp,zpsoil,raing,rainc,prcrate,                   &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 var3du,var3dv,var3dt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write ARPS history file using NetCDF 3.0 API.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2004/08/02
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'                                  ! Grid parameters
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
! Variable declaration
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: netid
  INTEGER, INTENT(IN)  :: itime                ! Time level, default 1.
  INTEGER, INTENT(IN)  :: packed               ! No pack implemented still, 0
  INTEGER, INTENT(IN)  :: nx,ny,nz,nzsoil,nstyps
  INTEGER, INTENT(IN)  :: grdbas               ! If this is a grid/base state array dump

  REAL,    INTENT(IN)  :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL,    INTENT(IN)  :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL,    INTENT(IN)  :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL,    INTENT(IN)  :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL,    INTENT(IN)  :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL,    INTENT(IN)  :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL,    INTENT(IN)  :: qscalar (nx,ny,nz,nscalar)

  REAL,    INTENT(IN)  :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL,    INTENT(IN)  :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                                               ! momentum. ( m**2/s )
  REAL,    INTENT(IN)  :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                                               ! momentum. ( m**2/s )
  REAL,    INTENT(IN)  :: ubar  (nx,ny,nz)     ! Base state x-velocity (m/s)
  REAL,    INTENT(IN)  :: vbar  (nx,ny,nz)     ! Base state y-velocity (m/s)
  REAL,    INTENT(IN)  :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL,    INTENT(IN)  :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL,    INTENT(IN)  :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL,    INTENT(IN)  :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                                               ! (kg/kg)
  REAL,    INTENT(IN)  :: x     (nx)           ! The x-coord. of the physical and
                                               ! computational grid. Defined at u-point.
  REAL,    INTENT(IN)  :: y     (ny)           ! The y-coord. of the physical and
                                               ! computational grid. Defined at v-point.
  REAL,    INTENT(IN)  :: z     (nz)           ! The z-coord. of the computational grid.
                                               ! Defined at w-point on the staggered grid.
  REAL,    INTENT(IN)  :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                                               ! w-point of the staggered grid.
  REAL,    INTENT(IN)  :: zpsoil (nx,ny,nzsoil)! The physical height coordinate defined at
                                               ! w-point of the soil.
  REAL,    INTENT(IN)  :: raing(nx,ny)         ! Grid supersaturation rain
  REAL,    INTENT(IN)  :: rainc(nx,ny)         ! Cumulus convective rain
  REAL,    INTENT(IN)  :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                                               ! prcrate(1,1,1) = total precip. rate
                                               ! prcrate(1,1,2) = grid scale precip. rate
                                               ! prcrate(1,1,3) = cumulus precip. rate
                                               ! prcrate(1,1,4) = microphysics precip. rate

  INTEGER, INTENT(IN)  :: soiltyp(nx,ny,nstyps)   ! Soil type
  REAL,    INTENT(IN)  :: stypfrct(nx,ny,nstyps)  ! Soil type fractions
  INTEGER, INTENT(IN)  :: vegtyp (nx,ny)          ! Vegetation type
  REAL,    INTENT(IN)  :: lai    (nx,ny)          ! Leaf Area Index
  REAL,    INTENT(IN)  :: roufns (nx,ny)          ! Surface roughness
  REAL,    INTENT(IN)  :: veg    (nx,ny)          ! Vegetation fraction

  REAL,    INTENT(IN)  :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL,    INTENT(IN)  :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL,    INTENT(IN)  :: wetcanp(nx,ny,0:nstyps)        ! Canopy water amount
  REAL,    INTENT(IN)  :: snowdpth(nx,ny)                ! Snow depth (m)

  REAL,    INTENT(IN)  :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL,    INTENT(IN)  :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL,    INTENT(IN)  :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL,    INTENT(IN)  :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL,    INTENT(IN)  :: radlwin(nx,ny)       ! Incominging longwave radiation

  REAL,    INTENT(IN)  :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL,    INTENT(IN)  :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL,    INTENT(IN)  :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL,    INTENT(IN)  :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL,  INTENT(INOUT) :: var3du(nx,  ny-1,nz-1)     ! Temporary work array
  REAL,  INTENT(INOUT) :: var3dv(nx-1,ny,  nz-1)     ! Temporary work array
  REAL,  INTENT(INOUT) :: var3dt(nx-1,ny-1,nz-1)     ! Temporary work array

!-----------------------------------------------------------------------
!
!  Local working arrays
!
!-----------------------------------------------------------------------

  REAL,    ALLOCATABLE :: var2d (:,:)
  REAL,    ALLOCATABLE :: var3dw(:,:,:)
  INTEGER, ALLOCATABLE :: var2di(:,:)
  INTEGER, ALLOCATABLE :: var3di(:,:,:)
  REAL,    ALLOCATABLE :: var4d (:,:,:,:)

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  CHARACTER(LEN=10), PARAMETER :: tmunit = 'seconds   '

  INTEGER :: zdim
  INTEGER :: i,j,k,is
  INTEGER :: nq
  INTEGER :: istatus

  INTEGER :: nxlg, nylg

  CHARACTER(LEN=4) :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  IF (mp_opt /= 1) joindmp = 0                ! Non-mpi run

  zdim = MAX(nz,nzsoil,nstyps+1,4)

  ALLOCATE(var2d (nx-1,ny-1),                 STAT = istatus)
  ALLOCATE(var2di(nx-1,ny-1),                 STAT = istatus)
  ALLOCATE(var3dw(nx-1,ny-1,zdim),            STAT = istatus)
  ALLOCATE(var3di(nx-1,ny-1,nstyps),          STAT = istatus)
  ALLOCATE(var4d (nx-1,ny-1,nzsoil,nstyps+1), STAT = istatus)

!-----------------------------------------------------------------------
!
! Define dimension, write global attribute and define variables,
! Only needed for the first time level
!
!-----------------------------------------------------------------------

  IF (myproc == 0) WRITE(6,'(/1x,a/)')                                  &
      'Defining NetCDF dimensions, global attributes and variables. '

  IF( itime == 1 ) THEN

!-----------------------------------------------------------------------
!
!  Define ARPS dimension and variables
!
!-----------------------------------------------------------------------

    nxlg = (nx-3)*nproc_x + 3
    nylg = (ny-3)*nproc_y + 3

    CALL net_define_dimension(netid,grdbas,nx,ny,nz,nzsoil,nstyps)
    CALL net_define_variables(netid,packed,grdbas,tmunit,nxlg,nylg,istatus)
                  ! nxlg, nylg used only in netwrt_general_att
  END IF

!-----------------------------------------------------------------------
!
! Beginning of writing variables
!
!-----------------------------------------------------------------------

  IF (grdbas == 1) THEN
    IF(myproc ==0) WRITE(6,'(1x,/,a/)')                                 &
                          'Writing history grid & base data.'
  ELSE
    IF(myproc ==0) WRITE(6,'(1x,/,a,f13.3/)')                           &
                          'Writing history data at time=', curtim
  END IF

  !
  ! Grid variables
  !
  IF (itime < 2 .AND. (grdout == 1 .OR. grdbas == 1) ) THEN

    CALL netwrt1d(netid,packed,0,'x_stag',x,nx)
    CALL netwrt1d(netid,packed,0,'y_stag',y,ny)
    CALL netwrt1d(netid,packed,0,'z_stag',z,nz)

    DO k = 1,nz
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3dw(i,j,k) = zp(i,j,k)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,0,'ZP',var3dw,nx-1,ny-1,nz)

    DO k = 1,nzsoil
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3dw(i,j,k) = zpsoil(i,j,k)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,0,'ZPSOIL',var3dw,nx-1,ny-1,nzsoil)

  END IF

  !
  ! Base state variables
  !
  IF(itime < 2 .AND. (basout == 1 .OR. grdbas == 1) ) THEN

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1, nx
          var3du(i,j,k) = ubar(i,j,k)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,0,'UBAR',var3du,nx,ny-1,nz-1)

    DO k = 1,nz-1
      DO j = 1,ny
        DO i = 1,nx-1
          var3dv(i,j,k) = vbar(i,j,k)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,0,'VBAR',var3dv,nx-1,ny,nz-1)

    var3dw(:,:,:) = 0.0
    CALL netwrt3d(netid,packed,0,'WBAR',var3dw,nx-1,ny-1,nz)

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3dt(i,j,k) = ptbar(i,j,k)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,0,'PTBAR',var3dt,nx-1,ny-1,nz-1)

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3dt(i,j,k) = pbar(i,j,k)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,0,'PBAR',var3dt,nx-1,ny-1,nz-1)

    IF (mstout == 1) THEN

      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3dt(i,j,k) = qvbar(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,0,'QVBAR',var3dt,nx-1,ny-1,nz-1)

    END IF

    IF (landout == 1) THEN

      DO k = 1,nstyps
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3di(i,j,k) = soiltyp(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3di(netid,packed,0,'SOILTYP',var3di,nx-1,ny-1,nstyps)

      DO k = 1,nstyps
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3dw(i,j,k) = stypfrct(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,0,'STYPFRCT',var3dw,nx-1,ny-1,nstyps)

      DO j = 1,ny-1
        DO i = 1,nx-1
          var2di(i,j) = vegtyp(i,j)
        END DO
      END DO
      CALL netwrt2di(netid,packed,0,'VEGTYP',var2di,nx-1,ny-1)

      DO j = 1,ny-1
        DO i = 1,nx-1
          var2d(i,j) = lai(i,j)
        END DO
      END DO
      CALL netwrt2d(netid,packed,0,'LAI',var2d,nx-1,ny-1)

      DO j = 1,ny-1
        DO i = 1,nx-1
          var2d(i,j) = roufns(i,j)
        END DO
      END DO
      CALL netwrt2d(netid,packed,0,'ROUFNS',var2d,nx-1,ny-1)

      DO j = 1,ny-1
        DO i = 1,nx-1
          var2d(i,j) = veg(i,j)
        END DO
      END DO
      CALL netwrt2d(netid,packed,0,'VEG',var2d,nx-1,ny-1)

    END IF
  END IF

  IF ( grdbas == 1 ) GOTO 3333

  CALL netwrtTime(netid,itime,'Time',curtim)

  IF (varout == 1) THEN

    IF (totout == 0) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            var3du(i,j,k)=u(i,j,k)-ubar(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'UPRT',var3du,nx,ny-1,nz-1)

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            var3dv(i,j,k)=v(i,j,k)-vbar(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'VPRT',var3dv,nx-1,ny,nz-1)

      DO k = 1,nz
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3dw(i,j,k) = w(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'WPRT',var3dw,nx-1,ny-1,nz)

      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3dt(i,j,k) = ptprt(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'PTPRT',var3dt,nx-1,ny-1,nz-1)

      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3dt(i,j,k) = pprt(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'PPRT',var3dt,nx-1,ny-1,nz-1)

    ELSE

      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx
            var3du(i,j,k) = u(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'U',var3du,nx,ny-1,nz-1)

      DO k = 1,nz-1
        DO j = 1,ny
          DO i = 1,nx-1
            var3dv(i,j,k) = v(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'V',var3dv,nx-1,ny,nz-1)

      DO k = 1,nz
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3dw(i,j,k) = w(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'W',var3dw,nx-1,ny-1,nz)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            var3dt(i,j,k) = ptbar(i,j,k) + ptprt(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'PT',var3dt,nx-1,ny-1,nz-1)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            var3dt(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'P',var3dt,nx-1,ny-1,nz-1)

    END IF     ! totout
  END IF    ! varout

  IF (mstout == 1) THEN

    IF (totout == 0) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            var3dt(i,j,k)=qv(i,j,k)-qvbar(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'QVPRT',var3dt,nx-1,ny-1,nz-1)

    ELSE

      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3dt(i,j,k) = qv(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'QV',var3dt,nx-1,ny-1,nz-1)

    END IF

    DO nq = 1,nscalar
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3dt(i,j,k) = qscalar(i,j,k,nq)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,upcase(qnames(nq)),var3dt,nx-1,ny-1,nz-1)
    END DO

    IF (rainout == 1) THEN

      DO j = 1,ny-1
        DO i = 1,nx-1
          var2d(i,j) = raing(i,j)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'RAING',var2d,nx-1,ny-1)

      DO j = 1,ny-1
        DO i = 1,nx-1
          var2d(i,j) = rainc(i,j)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'RAINC',var2d,nx-1,ny-1)

    END IF

    IF (prcout == 1) THEN

      DO k = 1,4
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3dw(i,j,k) = prcrate(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'PRCRATE1',var3dw(:,:,1),nx-1,ny-1)
      CALL netwrt2d(netid,packed,itime,'PRCRATE2',var3dw(:,:,2),nx-1,ny-1)
      CALL netwrt2d(netid,packed,itime,'PRCRATE3',var3dw(:,:,3),nx-1,ny-1)
      CALL netwrt2d(netid,packed,itime,'PRCRATE4',var3dw(:,:,4),nx-1,ny-1)

    END IF

  END IF      ! mstout == 1

  IF (tkeout == 1) THEN

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3dt(i,j,k) = tke(i,j,k)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,itime,'TKE',var3dt,nx-1,ny-1,nz-1)

  END IF

  IF (trbout == 1) THEN

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3dt(i,j,k) = kmh(i,j,k)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,itime,'KMH',var3dt,nx-1,ny-1,nz-1)

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3dt(i,j,k) = kmv(i,j,k)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,itime,'KMV',var3dt,nx-1,ny-1,nz-1)

  END IF  ! trbout

  IF (sfcout == 1) THEN

    DO is = 0,nstyps
      DO k = 1,nzsoil
        DO j = 1,ny-1
          DO i = 1,nx-1
            var4d(i,j,k,is+1) = tsoil(i,j,k,is)
          END DO
        END DO
      END DO
    END DO
    CALL netwrt4d(netid,packed,itime,'TSOIL',var4d,nx-1,ny-1,nzsoil,nstyps+1)

    DO is = 0,nstyps
      DO k = 1,nzsoil
        DO j = 1,ny-1
          DO i = 1,nx-1
            var4d(i,j,k,is+1) = qsoil(i,j,k,is)
          END DO
        END DO
      END DO
    END DO
    CALL netwrt4d(netid,packed,itime,'QSOIL',var4d,nx-1,ny-1,nzsoil,nstyps+1)

    DO is = 0,nstyps
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3dw(i,j,is+1) = wetcanp(i,j,is)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,itime,'WETCANP',var3dw,nx-1,ny-1,nstyps+1)

    IF (snowout == 1) THEN

      DO j = 1,ny-1
        DO i = 1,nx-1
          var2d(i,j) = snowdpth(i,j)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'SNOWDPTH',var2d,nx-1,ny-1)

    END IF

  END IF  ! sfcout

  IF (radout == 1) THEN

    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3dt(i,j,k) = radfrc(i,j,k)
        END DO
      END DO
    END DO
    CALL netwrt3d(netid,packed,itime,'RADFRC',var3dt,nx-1,ny-1,nz-1)

    DO j = 1,ny-1
      DO i = 1,nx-1
        var2d(i,j) = radsw(i,j)
      END DO
    END DO
    CALL netwrt2d(netid,packed,itime,'RADSW',var2d,nx-1,ny-1)

    DO j = 1,ny-1
      DO i = 1,nx-1
        var2d(i,j) = rnflx(i,j)
      END DO
    END DO
    CALL netwrt2d(netid,packed,itime,'RNFLX',var2d,nx-1,ny-1)

    DO j = 1,ny-1
      DO i = 1,nx-1
        var2d(i,j) = radswnet(i,j)
      END DO
    END DO
    CALL netwrt2d(netid,packed,itime,'RADSWNET',var2d,nx-1,ny-1)

    DO j = 1,ny-1
      DO i = 1,nx-1
        var2d(i,j) = radlwin(i,j)
      END DO
    END DO
    CALL netwrt2d(netid,packed,itime,'RADLWIN',var2d,nx-1,ny-1)

  END IF

  IF (flxout == 1) THEN

    DO j = 1,ny-1
      DO i = 1,nx-1
        var2d(i,j) = usflx(i,j)
      END DO
    END DO
    CALL netwrt2d(netid,packed,itime,'USFLX',var2d,nx-1,ny-1)

    DO j = 1,ny-1
      DO i = 1,nx-1
        var2d(i,j) = vsflx(i,j)
      END DO
    END DO
    CALL netwrt2d(netid,packed,itime,'VSFLX',var2d,nx-1,ny-1)

    DO j = 1,ny-1
      DO i = 1,nx-1
        var2d(i,j) = ptsflx(i,j)
      END DO
    END DO
    CALL netwrt2d(netid,packed,itime,'PTSFLX',var2d,nx-1,ny-1)

    DO j = 1,ny-1
      DO i = 1,nx-1
        var2d(i,j) = qvsflx(i,j)
      END DO
    END DO
    CALL netwrt2d(netid,packed,itime,'QVSFLX',var2d,nx-1,ny-1)

  END IF

  3333 CONTINUE

  DEALLOCATE(var2d, var2di)
  DEALLOCATE(var3dw,var3di)
  DEALLOCATE(var4d)

  RETURN
END SUBROUTINE netdump
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE NETJOINDUMP                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE netjoindump(netid,itime,packed,nx,ny,nz,nzsoil,nstyps,grdbas,    &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 kmh,kmv,ubar,vbar,ptbar,pbar,rhobar,qvbar,             &
                 x,y,z,zp,zpsoil,raing,rainc,prcrate,                   &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write ARPS history file using NetCDF 3.0 API.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2004/08/02
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'                                  ! Grid parameters
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
! Variable declaration
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: netid
  INTEGER, INTENT(IN)  :: itime                ! Time level, default 1.
  INTEGER, INTENT(IN)  :: packed               ! No pack implemented still, 0
  INTEGER, INTENT(IN)  :: nx,ny,nz,nzsoil,nstyps
  INTEGER, INTENT(IN)  :: grdbas               ! If this is a grid/base state array dump

  REAL,    INTENT(IN)  :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL,    INTENT(IN)  :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL,    INTENT(IN)  :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL,    INTENT(IN)  :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL,    INTENT(IN)  :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL,    INTENT(IN)  :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL,    INTENT(IN)  :: qscalar  (nx,ny,nz,nscalar)

  REAL,    INTENT(IN)  :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL,    INTENT(IN)  :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                                               ! momentum. ( m**2/s )
  REAL,    INTENT(IN)  :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                                               ! momentum. ( m**2/s )
  REAL,    INTENT(IN)  :: ubar  (nx,ny,nz)     ! Base state x-velocity (m/s)
  REAL,    INTENT(IN)  :: vbar  (nx,ny,nz)     ! Base state y-velocity (m/s)
  REAL,    INTENT(IN)  :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL,    INTENT(IN)  :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL,    INTENT(IN)  :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL,    INTENT(IN)  :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                                               ! (kg/kg)
  REAL,    INTENT(IN)  :: x     (nx)           ! The x-coord. of the physical and
                                               ! computational grid. Defined at u-point.
  REAL,    INTENT(IN)  :: y     (ny)           ! The y-coord. of the physical and
                                               ! computational grid. Defined at v-point.
  REAL,    INTENT(IN)  :: z     (nz)           ! The z-coord. of the computational grid.
                                               ! Defined at w-point on the staggered grid.
  REAL,    INTENT(IN)  :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                                               ! w-point of the staggered grid.
  REAL,    INTENT(IN)  :: zpsoil (nx,ny,nzsoil)! The physical height coordinate defined at
                                               ! w-point of the soil.
  REAL,    INTENT(IN)  :: raing(nx,ny)         ! Grid supersaturation rain
  REAL,    INTENT(IN)  :: rainc(nx,ny)         ! Cumulus convective rain
  REAL,    INTENT(IN)  :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                                               ! prcrate(1,1,1) = total precip. rate
                                               ! prcrate(1,1,2) = grid scale precip. rate
                                               ! prcrate(1,1,3) = cumulus precip. rate
                                               ! prcrate(1,1,4) = microphysics precip. rate

  INTEGER, INTENT(IN)  :: soiltyp(nx,ny,nstyps)   ! Soil type
  REAL,    INTENT(IN)  :: stypfrct(nx,ny,nstyps)  ! Soil type fractions
  INTEGER, INTENT(IN)  :: vegtyp (nx,ny)          ! Vegetation type
  REAL,    INTENT(IN)  :: lai    (nx,ny)          ! Leaf Area Index
  REAL,    INTENT(IN)  :: roufns (nx,ny)          ! Surface roughness
  REAL,    INTENT(IN)  :: veg    (nx,ny)          ! Vegetation fraction

  REAL,    INTENT(IN)  :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL,    INTENT(IN)  :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL,    INTENT(IN)  :: wetcanp(nx,ny,0:nstyps)        ! Canopy water amount
  REAL,    INTENT(IN)  :: snowdpth(nx,ny)                ! Snow depth (m)

  REAL,    INTENT(IN)  :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL,    INTENT(IN)  :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL,    INTENT(IN)  :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL,    INTENT(IN)  :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL,    INTENT(IN)  :: radlwin(nx,ny)       ! Incominging longwave radiation

  REAL,    INTENT(IN)  :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL,    INTENT(IN)  :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL,    INTENT(IN)  :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL,    INTENT(IN)  :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL,  INTENT(INOUT) :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL,  INTENT(INOUT) :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL,  INTENT(INOUT) :: tem3  (nx,ny,nz)     ! Temporary work array

!-----------------------------------------------------------------------
!
!  Local working arrays
!
!-----------------------------------------------------------------------

  REAL,    ALLOCATABLE :: out1d(:)
  REAL,    ALLOCATABLE :: out3d(:,:,:)
  REAL,    ALLOCATABLE :: out4d(:,:,:,:)
  INTEGER, ALLOCATABLE :: out3di(:,:,:)

  REAL,    ALLOCATABLE :: var3du(:,:,:)
  REAL,    ALLOCATABLE :: var3dv(:,:,:)
  REAL,    ALLOCATABLE :: var3dw(:,:,:)
  REAL,    ALLOCATABLE :: var4d (:,:,:,:)
  INTEGER, ALLOCATABLE :: var3di(:,:,:)

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  CHARACTER(LEN=10), PARAMETER :: tmunit = 'seconds   '

  INTEGER :: nxlg,nylg
  INTEGER :: i,j,k,is
  INTEGER :: nq
  INTEGER :: istatus

  CHARACTER(LEN=4) :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

  ALLOCATE(out1d (MAX(nxlg,nylg)),                    STAT = istatus)
  ALLOCATE(out3d (nxlg,nylg,MAX(nz,nzsoil,nstyps+1)), STAT = istatus)
  ALLOCATE(out3di(nxlg,nylg,nstyps),                  STAT = istatus)
  ALLOCATE(out4d (nxlg,nylg,nzsoil,nstyps+1),         STAT = istatus)

  IF (myproc == 0) THEN
    ALLOCATE(var3du(nxlg,  nylg-1,nz-1),                    STAT = istatus)
    ALLOCATE(var3dv(nxlg-1,nylg,  nz-1),                    STAT = istatus)
    ALLOCATE(var3dw(nxlg-1,nylg-1,MAX(nz,nstyps+1,nzsoil)), STAT = istatus)
    ALLOCATE(var3di(nxlg-1,nylg-1,nstyps),                  STAT = istatus)
    ALLOCATE(var4d (nxlg-1,nylg-1,nzsoil,nstyps+1),         STAT = istatus)
  END IF

!-----------------------------------------------------------------------
!
! Define dimension, write global attribute and define variables,
! Only needed for the first time level
!
!-----------------------------------------------------------------------

  IF (myproc == 0) WRITE(6,'(/1x,a/)')                                  &
      'Defining NetCDF dimensions, global attribute and variables. '

  IF( (itime == 1) .AND. (myproc == 0) ) THEN
!
!-----------------------------------------------------------------------
!
!  Define ARPS dimension and variables
!
!-----------------------------------------------------------------------

    CALL net_define_dimension(netid,grdbas,nxlg,nylg,nz,nzsoil,nstyps)
    CALL net_define_variables(netid,packed,grdbas,tmunit,nxlg,nylg,istatus)
  END IF

!-----------------------------------------------------------------------
!
! Beginning of writing variables
!
!-----------------------------------------------------------------------

  IF (grdbas == 1) THEN
    IF(myproc ==0) WRITE(6,'(1x,/,a/)')                                 &
                          'Writing history grid & base data.'
  ELSE
    IF(myproc ==0) WRITE(6,'(1x,/,a,f13.3/)')                           &
                          'Writing history data at time=', curtim
  END IF

  !
  ! Grid variables
  !
  IF (itime == 1 .AND. (grdout == 1 .OR. grdbas == 1) ) THEN

    CALL mpimerge1dx(x,nx,out1d)
    IF( myproc == 0 ) THEN
      CALL netwrt1d(netid,packed,0,'x_stag',out1d,nxlg)
    END IF

    CALL mpimerge1dy(y,ny,out1d)
    IF( myproc == 0  ) THEN
      CALL netwrt1d(netid,packed,0,'y_stag',out1d,nylg)
    END IF

    IF( myproc == 0  ) THEN
      CALL netwrt1d(netid,packed,0,'z_stag',z,nz)
    END IF

    CALL mpimerge3d(zp,nx,ny,nz,out3d)
    IF( myproc == 0  ) THEN
      DO k = 1,nz
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,k) = out3d(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,0,'ZP',var3dw,nxlg-1,nylg-1,nz)
    END IF

    CALL mpimerge3d(zpsoil,nx,ny,nzsoil,out3d)
    IF( myproc == 0  ) THEN
      DO k = 1,nzsoil
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,k) = out3d(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,0,'ZPSOIL',var3dw,nxlg-1,nylg-1,nzsoil)
    END IF

  END IF

  !
  ! Base state variables
  !
  IF(itime == 1 .AND. (basout == 1 .OR. grdbas == 1) ) THEN

    CALL mpimerge3d(ubar,nx,ny,nz,out3d)
    IF( myproc == 0  ) THEN
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg
            var3du(i,j,k) = out3d(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,0,'UBAR',var3du,nxlg,nylg-1,nz-1)
    END IF

    CALL mpimerge3d(vbar,nx,ny,nz,out3d)
    IF( myproc == 0  ) THEN
      DO k = 1,nz-1
        DO j = 1,nylg
          DO i = 1,nxlg-1
            var3dv(i,j,k) = out3d(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,0,'VBAR',var3dv,nxlg-1,nylg,nz-1)
    END IF

    IF( myproc == 0  ) THEN
      DO k = 1,nz
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,k) = 0.0
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,0,'WBAR',var3dw,nxlg-1,nylg-1,nz)
    END IF

    CALL mpimerge3d(ptbar,nx,ny,nz,out3d)
    IF( myproc == 0  ) THEN
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,k) = out3d(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,0,'PTBAR',var3dw,nxlg-1,nylg-1,nz-1)
    END IF

    CALL mpimerge3d(pbar,nx,ny,nz,out3d)
    IF( myproc == 0  ) THEN
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,k) = out3d(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,0,'PBAR',var3dw,nxlg-1,nylg-1,nz-1)
    END IF

    IF (mstout == 1) THEN

      CALL mpimerge3d(qvbar,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,0,'QVBAR',var3dw,nxlg-1,nylg-1,nz-1)
      END IF

    END IF

    IF (landout == 1) THEN

      CALL mpimerge3di(soiltyp,nx,ny,nstyps,out3di)
      IF( myproc == 0  ) THEN
        DO is = 1,nstyps
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3di(i,j,is) = out3di(i,j,is)
            END DO
          END DO
        END DO
        CALL netwrt3di(netid,packed,0,'SOILTYP',var3di,nxlg-1,nylg-1,nstyps)
      END IF

      CALL mpimerge3d(stypfrct,nx,ny,nstyps,out3d)
      IF( myproc == 0  ) THEN
        DO is = 1,nstyps
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,is) = out3d(i,j,is)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,0,'STYPFRCT',var3dw,nxlg-1,nylg-1,nstyps)
      END IF

      CALL mpimerge2di(vegtyp,nx,ny,out3di)
      IF( myproc == 0  ) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3di(i,j,1) = out3di(i,j,1)
          END DO
        END DO
        CALL netwrt2di(netid,packed,0,'VEGTYP',var3di,nxlg-1,nylg-1)
      END IF

      CALL mpimerge2d(lai,nx,ny,out3d)
      IF( myproc == 0  ) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,1) = out3d(i,j,1)
          END DO
        END DO
        CALL netwrt2d(netid,packed,0,'LAI',var3dw,nxlg-1,nylg-1)
      END IF

      CALL mpimerge2d(roufns,nx,ny,out3d)
      IF( myproc == 0  ) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,1) = out3d(i,j,1)
          END DO
        END DO
        CALL netwrt2d(netid,packed,0,'ROUFNS',var3dw,nxlg-1,nylg-1)
      END IF

      CALL mpimerge2d(veg,nx,ny,out3d)
      IF( myproc == 0  ) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,1) = out3d(i,j,1)
          END DO
        END DO
        CALL netwrt2d(netid,packed,0,'VEG',var3dw,nxlg-1,nylg-1)
      END IF

    END IF
  END IF

  IF ( grdbas == 1 ) GOTO 3333

  IF( myproc == 0  ) CALL netwrtTime(netid,itime,'Time',curtim)

  IF (varout == 1) THEN

    IF (totout == 0) THEN

      tem1(:,:,:) = u(:,:,:) - ubar(:,:,:)
      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg
              var3du(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'UPRT',var3du,nxlg,nylg-1,nz-1)
      END IF

      tem1(:,:,:) = v(:,:,:) - vbar(:,:,:)
      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg
            DO i = 1,nxlg-1
              var3dv(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'VPRT',var3dv,nxlg-1,nylg,nz-1)
      END IF

      CALL mpimerge3d(w,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'WPRT',var3dw,nxlg-1,nylg-1,nz)
      END IF

      CALL mpimerge3d(ptprt,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'PTPRT',var3dw,nxlg-1,nylg-1,nz-1)
      END IF

      CALL mpimerge3d(pprt,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'PPRT',var3dw,nxlg-1,nylg-1,nz-1)
      END IF

    ELSE

      CALL mpimerge3d(u,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg
              var3du(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'U',var3du,nxlg,nylg-1,nz-1)
      END IF

      CALL mpimerge3d(v,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg
            DO i = 1,nxlg-1
              var3dv(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'V',var3dv,nxlg-1,nylg,nz-1)
      END IF

      CALL mpimerge3d(w,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'W',var3dw,nxlg-1,nylg-1,nz)
      END IF

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = ptbar(i,j,k) + ptprt(i,j,k)
          END DO
        END DO
      END DO

      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'PT',var3dw,nxlg-1,nylg-1,nz-1)
      END IF

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
          END DO
        END DO
      END DO

      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'P',var3dw,nxlg-1,nylg-1,nz-1)
      END IF

    END IF     ! totout
  END IF    ! varout

  IF (mstout == 1) THEN

    IF (totout == 0) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=qv(i,j,k)-qvbar(i,j,k)
          END DO
        END DO
      END DO

      CALL mpimerge3d(tem1,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'QVPRT',var3dw,nxlg-1,nylg-1,nz-1)
      END IF

    ELSE

      CALL mpimerge3d(qv,nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'QV',var3dw,nxlg-1,nylg-1,nz-1)
      END IF

    END IF

    DO nq = 1,nscalar
      CALL mpimerge3d(qscalar(:,:,:,nq),nx,ny,nz,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,nz-1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,upcase(qnames(nq)),var3dw,nxlg-1,nylg-1,nz-1)
      END IF
    END DO

    IF (rainout == 1) THEN

      CALL mpimerge2d(raing,nx,ny,out3d)
      IF( myproc == 0  ) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,1) = out3d(i,j,1)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'RAING',var3dw,nxlg-1,nylg-1)
      END IF

      CALL mpimerge2d(rainc,nx,ny,out3d)
      IF( myproc == 0  ) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,1) = out3d(i,j,1)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'RAINC',var3dw,nxlg-1,nylg-1)
      END IF

    END IF

    IF (prcout == 1) THEN

      CALL mpimerge3d(prcrate,nx,ny,4,out3d)
      IF( myproc == 0  ) THEN
        DO k = 1,4
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3dw(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'PRCRATE1',var3dw(:,:,1),nxlg-1,nylg-1)
        CALL netwrt2d(netid,packed,itime,'PRCRATE2',var3dw(:,:,2),nxlg-1,nylg-1)
        CALL netwrt2d(netid,packed,itime,'PRCRATE3',var3dw(:,:,3),nxlg-1,nylg-1)
        CALL netwrt2d(netid,packed,itime,'PRCRATE4',var3dw(:,:,4),nxlg-1,nylg-1)
      END IF

    END IF

  END IF      ! mstout == 1

  IF (tkeout == 1) THEN

    CALL mpimerge3d(tke,nx,ny,nz,out3d)
    IF( myproc == 0  ) THEN
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,k) = out3d(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'TKE',var3dw,nxlg-1,nylg-1,nz-1)
    END IF

  END IF

  IF (trbout == 1) THEN

    CALL mpimerge3d(kmh,nx,ny,nz,out3d)
    IF( myproc == 0  ) THEN
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,k) = out3d(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'KMH',var3dw,nxlg-1,nylg-1,nz-1)
    END IF

    CALL mpimerge3d(kmv,nx,ny,nz,out3d)
    IF( myproc == 0  ) THEN
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,k) = out3d(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'KMV',var3dw,nxlg-1,nylg-1,nz-1)
    END IF

  END IF  ! trbout

  IF (sfcout == 1) THEN

    CALL mpimerge4d(tsoil,nx,ny,nzsoil,nstyps+1,out4d)
    IF( myproc == 0 ) THEN
      DO is = 1, nstyps+1
        DO k = 1, nzsoil
          DO j = 1, nylg-1
            DO i = 1, nxlg-1
              var4d(i,j,k,is) = out4d(i,j,k,is)
            END DO
          END DO
        END DO
      END DO
      CALL netwrt4d(netid,packed,itime,'TSOIL',var4d,nxlg-1,nylg-1,nzsoil,nstyps+1)
    END IF

    CALL mpimerge4d(qsoil,nx,ny,nzsoil,nstyps+1,out4d)
    IF( myproc == 0  ) THEN
      DO is = 1, nstyps+1
        DO k = 1, nzsoil
          DO j = 1, nylg-1
            DO i = 1, nxlg-1
              var4d(i,j,k,is) = out4d(i,j,k,is)
            END DO
          END DO
        END DO
      END DO
      CALL netwrt4d(netid,packed,itime,'QSOIL',var4d,nxlg-1,nylg-1,nzsoil,nstyps+1)
    END IF

    CALL mpimerge3d(wetcanp,nx,ny,nstyps+1,out3d)
    IF( myproc == 0  ) THEN
      DO is = 1, nstyps+1
        DO j = 1, nylg-1
          DO i = 1, nxlg-1
            var3dw(i,j,is) = out3d(i,j,is)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'WETCANP',var3dw,nxlg-1,nylg-1,nstyps+1)
    END IF

    IF (snowout == 1) THEN

      CALL mpimerge2d(snowdpth,nx,ny,out3d)
      IF( myproc == 0  ) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,1) = out3d(i,j,1)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'SNOWDPTH',var3dw,nxlg-1,nylg-1)
      END IF

    END IF

  END IF  ! sfcout

  IF (radout == 1) THEN

    CALL mpimerge3d(radfrc,nx,ny,nz,out3d)
    IF( myproc == 0  ) THEN
      DO k = 1,nz-1
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3dw(i,j,k) = out3d(i,j,k)
          END DO
        END DO
      END DO
      CALL netwrt3d(netid,packed,itime,'RADFRC',var3dw,nxlg-1,nylg-1,nz-1)
    END IF

    CALL mpimerge2d(radsw,nx,ny,out3d)
    IF( myproc == 0  ) THEN
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3dw(i,j,1) = out3d(i,j,1)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'RADSW',var3dw,nxlg-1,nylg-1)
    END IF

    CALL mpimerge2d(rnflx,nx,ny,out3d)
    IF( myproc == 0  ) THEN
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3dw(i,j,1) = out3d(i,j,1)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'RNFLX',var3dw,nxlg-1,nylg-1)
    END IF

    CALL mpimerge2d(radswnet,nx,ny,out3d)
    IF( myproc == 0  ) THEN
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3dw(i,j,1) = out3d(i,j,1)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'RADSWNET',var3dw,nxlg-1,nylg-1)
    END IF

    CALL mpimerge2d(radlwin,nx,ny,out3d)
    IF( myproc == 0  ) THEN
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3dw(i,j,1) = out3d(i,j,1)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'RADLWIN',var3dw,nxlg-1,nylg-1)
    END IF

  END IF

  IF (flxout == 1) THEN

    CALL mpimerge2d(usflx,nx,ny,out3d)
    IF( myproc == 0  ) THEN
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3dw(i,j,1) = out3d(i,j,1)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'USFLX',var3dw,nxlg-1,nylg-1)
    END IF

    CALL mpimerge2d(vsflx,nx,ny,out3d)
    IF( myproc == 0  ) THEN
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3dw(i,j,1) = out3d(i,j,1)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'VSFLX',var3dw,nxlg-1,nylg-1)
    END IF

    CALL mpimerge2d(ptsflx,nx,ny,out3d)
    IF( myproc == 0  ) THEN
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3dw(i,j,1) = out3d(i,j,1)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'PTSFLX',var3dw,nxlg-1,nylg-1)
    END IF

    CALL mpimerge2d(qvsflx,nx,ny,out3d)
    IF( myproc == 0  ) THEN
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3dw(i,j,1) = out3d(i,j,1)
        END DO
      END DO
      CALL netwrt2d(netid,packed,itime,'QVSFLX',var3dw,nxlg-1,nylg-1)
    END IF

  END IF

  3333 CONTINUE

  DEALLOCATE(out1d,out3d,out4d,out3di)
  IF (myproc == 0) DEALLOCATE(var3du,var3dv,var3dw,var4d,var3di)

  RETURN
END SUBROUTINE netjoindump
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE NETSPLITDUMP               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE netsplitdump(filename,itime,packed,nx,ny,nz,nzsoil,nstyps,grdbas,    &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 kmh,kmv,ubar,vbar,ptbar,pbar,rhobar,qvbar,             &
                 x,y,z,zp,zpsoil,raing,rainc,prcrate,                   &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write ARPS history file using NetCDF 3.0 API.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2004/08/02
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'                                  ! Grid parameters
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------

  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER, INTENT(IN)  :: itime                ! Time level, default 1.
  INTEGER, INTENT(IN)  :: packed               ! No pack implemented still, 0
  INTEGER, INTENT(IN)  :: nx,ny,nz,nzsoil,nstyps
  INTEGER, INTENT(IN)  :: grdbas               ! If this is a grid/base state array dump

  REAL,    INTENT(IN)  :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL,    INTENT(IN)  :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL,    INTENT(IN)  :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL,    INTENT(IN)  :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL,    INTENT(IN)  :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL,    INTENT(IN)  :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL,    INTENT(IN)  :: qscalar    (nx,ny,nz,nscalar)

  REAL,    INTENT(IN)  :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL,    INTENT(IN)  :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                                               ! momentum. ( m**2/s )
  REAL,    INTENT(IN)  :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                                               ! momentum. ( m**2/s )
  REAL,    INTENT(IN)  :: ubar  (nx,ny,nz)     ! Base state x-velocity (m/s)
  REAL,    INTENT(IN)  :: vbar  (nx,ny,nz)     ! Base state y-velocity (m/s)
  REAL,    INTENT(IN)  :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL,    INTENT(IN)  :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL,    INTENT(IN)  :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL,    INTENT(IN)  :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                                               ! (kg/kg)
  REAL,    INTENT(IN)  :: x     (nx)           ! The x-coord. of the physical and
                                               ! computational grid. Defined at u-point.
  REAL,    INTENT(IN)  :: y     (ny)           ! The y-coord. of the physical and
                                               ! computational grid. Defined at v-point.
  REAL,    INTENT(IN)  :: z     (nz)           ! The z-coord. of the computational grid.
                                               ! Defined at w-point on the staggered grid.
  REAL,    INTENT(IN)  :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                                               ! w-point of the staggered grid.
  REAL,    INTENT(IN)  :: zpsoil (nx,ny,nzsoil)! The physical height coordinate defined at
                                               ! w-point of the soil.
  REAL,    INTENT(IN)  :: raing(nx,ny)         ! Grid supersaturation rain
  REAL,    INTENT(IN)  :: rainc(nx,ny)         ! Cumulus convective rain
  REAL,    INTENT(IN)  :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                                               ! prcrate(1,1,1) = total precip. rate
                                               ! prcrate(1,1,2) = grid scale precip. rate
                                               ! prcrate(1,1,3) = cumulus precip. rate
                                               ! prcrate(1,1,4) = microphysics precip. rate

  INTEGER, INTENT(IN)  :: soiltyp(nx,ny,nstyps)   ! Soil type
  REAL,    INTENT(IN)  :: stypfrct(nx,ny,nstyps)  ! Soil type fractions
  INTEGER, INTENT(IN)  :: vegtyp (nx,ny)          ! Vegetation type
  REAL,    INTENT(IN)  :: lai    (nx,ny)          ! Leaf Area Index
  REAL,    INTENT(IN)  :: roufns (nx,ny)          ! Surface roughness
  REAL,    INTENT(IN)  :: veg    (nx,ny)          ! Vegetation fraction

  REAL,    INTENT(IN)  :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL,    INTENT(IN)  :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL,    INTENT(IN)  :: wetcanp(nx,ny,0:nstyps)        ! Canopy water amount
  REAL,    INTENT(IN)  :: snowdpth(nx,ny)                ! Snow depth (m)

  REAL,    INTENT(IN)  :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL,    INTENT(IN)  :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL,    INTENT(IN)  :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL,    INTENT(IN)  :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL,    INTENT(IN)  :: radlwin(nx,ny)       ! Incominging longwave radiation

  REAL,    INTENT(IN)  :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL,    INTENT(IN)  :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL,    INTENT(IN)  :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL,    INTENT(IN)  :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  !REAL,  INTENT(INOUT) :: var3du(nx,  ny-1,nz-1)     ! Temporary work array
  !REAL,  INTENT(INOUT) :: var3dv(nx-1,ny,  nz-1)     ! Temporary work array
  !REAL,  INTENT(INOUT) :: var3dt(nx-1,ny-1,nz-1)     ! Temporary work array

!-----------------------------------------------------------------------
!
!  Local working arrays
!
!-----------------------------------------------------------------------

  REAL,    ALLOCATABLE :: var2d (:,:)
  REAL,    ALLOCATABLE :: var3du(:,:,:),var3dv(:,:,:),var3dt(:,:,:),var3dw(:,:,:)
  INTEGER, ALLOCATABLE :: var2di(:,:)
  INTEGER, ALLOCATABLE :: var3di(:,:,:)
  REAL,    ALLOCATABLE :: var4d (:,:,:,:)

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  CHARACTER(LEN=10), PARAMETER :: tmunit = 'seconds   '

  INTEGER :: zdim
  INTEGER :: i,j,k,is
  INTEGER :: nq
  INTEGER :: istatus

  INTEGER :: nxlg, nylg

  INTEGER :: netid
  INTEGER :: npxout, npyout
  INTEGER :: nxout,  nyout
  INTEGER :: ipx,    jpy
  INTEGER :: ia,     ja

  INTEGER :: lenbase
  CHARACTER(LEN=256) :: outfilename

  CHARACTER(LEN=4) :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  IF (mp_opt /= 1) joindmp = 0                ! Non-mpi run

  npxout = nproc_x_out/nproc_x
  npyout = nproc_y_out/nproc_y

  nxout = (nx-3)/npxout + 3
  nyout = (ny-3)/npyout + 3

  lenbase = INDEX(filename,'_',.TRUE.)
  IF (lenbase < 0) lenbase = LEN_TRIM(filename)

  zdim = MAX(nz,nzsoil,nstyps+1,4)

  ALLOCATE(var2d (nxout-1,nyout-1),                 STAT = istatus)
  ALLOCATE(var2di(nxout-1,nyout-1),                 STAT = istatus)
  ALLOCATE(var3du(nxout,  nyout-1,nz-1),            STAT = istatus)
  ALLOCATE(var3dv(nxout-1,nyout,  nz-1),            STAT = istatus)
  ALLOCATE(var3dw(nxout-1,nyout-1,zdim),            STAT = istatus)
  ALLOCATE(var3dt(nxout-1,nyout-1,nz-1),            STAT = istatus)
  ALLOCATE(var3di(nxout-1,nyout-1,nstyps),          STAT = istatus)
  ALLOCATE(var4d (nxout-1,nyout-1,nzsoil,nstyps+1), STAT = istatus)

!-----------------------------------------------------------------------
!
! Define dimension, write global attribute and define variables,
! Only needed for the first time level
!
!-----------------------------------------------------------------------

  IF (myproc == 0) WRITE(6,'(/1x,a/)')                                  &
      'Defining NetCDF dimensions, global attribute and variables. '

  DO jpy = 1, npyout
    DO ipx = 1, npxout

      ia = (ipx-1)*(nxout-3)
      ja = (jpy-1)*(nyout-3)

      CALL gtsplitfn(filename(1:lenbase-1),npxout,npyout,loc_x,loc_y,ipx,jpy,  &
                     0,0,0,lvldbg,outfilename,istatus)

      CALL netopen(TRIM(outfilename),'C',netid)

!-----------------------------------------------------------------------
!
!  Define ARPS dimension and variables
!
!-----------------------------------------------------------------------

      nxlg = (nx-3)*nproc_x + 3
      nylg = (ny-3)*nproc_y + 3

      CALL net_define_dimension(netid,grdbas,nxout,nyout,nz,nzsoil,nstyps)
      CALL net_define_variables(netid,packed,grdbas,tmunit,nxlg,nylg,istatus)
                    ! nxlg, nylg used only in netwrt_general_att

!-----------------------------------------------------------------------
!
! Beginning of writing variables
!
!-----------------------------------------------------------------------

      IF (grdbas == 1) THEN
        IF(myproc ==0) WRITE(6,'(1x,/,a/)')                                 &
                              'Writing history grid & base data.'
      ELSE
        IF(myproc ==0) WRITE(6,'(1x,/,a,f13.3/)')                           &
                              'Writing history data at time=', curtim
      END IF

      !
      ! Grid variables
      !
      IF (itime < 2 .AND. (grdout == 1 .OR. grdbas == 1) ) THEN

        CALL netwrt1d(netid,packed,0,'x_stag',x(ia+1),nxout)
        CALL netwrt1d(netid,packed,0,'y_stag',y(ja+1),nyout)
        CALL netwrt1d(netid,packed,0,'z_stag',z,nz)

        DO k = 1,nz
          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var3dw(i,j,k) = zp(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,0,'ZP',var3dw,nxout-1,nyout-1,nz)

        DO k = 1,nzsoil
          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var3dw(i,j,k) = zpsoil(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,0,'ZPSOIL',var3dw,nxout-1,nyout-1,nzsoil)

      END IF

      !
      ! Base state variables
      !
      IF(itime < 2 .AND. (basout == 1 .OR. grdbas == 1) ) THEN

        DO k = 1,nz-1
          DO j = 1,nyout-1
            DO i = 1, nxout
              var3du(i,j,k) = ubar(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,0,'UBAR',var3du,nxout,nyout-1,nz-1)

        DO k = 1,nz-1
          DO j = 1,nyout
            DO i = 1,nxout-1
              var3dv(i,j,k) = vbar(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,0,'VBAR',var3dv,nxout-1,nyout,nz-1)

        var3dw(:,:,:) = 0.0
        CALL netwrt3d(netid,packed,0,'WBAR',var3dw,nxout-1,nyout-1,nz)

        DO k = 1,nz-1
          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var3dt(i,j,k) = ptbar(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,0,'PTBAR',var3dt,nxout-1,nyout-1,nz-1)

        DO k = 1,nz-1
          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var3dt(i,j,k) = pbar(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,0,'PBAR',var3dt,nxout-1,nyout-1,nz-1)

        IF (mstout == 1) THEN

          DO k = 1,nz-1
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dt(i,j,k) = qvbar(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,0,'QVBAR',var3dt,nxout-1,nyout-1,nz-1)

        END IF

        IF (landout == 1) THEN

          DO k = 1,nstyps
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3di(i,j,k) = soiltyp(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3di(netid,packed,0,'SOILTYP',var3di,nxout-1,nyout-1,nstyps)

          DO k = 1,nstyps
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dw(i,j,k) = stypfrct(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,0,'STYPFRCT',var3dw,nxout-1,nyout-1,nstyps)

          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var2di(i,j) = vegtyp(i+ia,j+ja)
            END DO
          END DO
          CALL netwrt2di(netid,packed,0,'VEGTYP',var2di,nxout-1,nyout-1)

          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var2d(i,j) = lai(i+ia,j+ja)
            END DO
          END DO
          CALL netwrt2d(netid,packed,0,'LAI',var2d,nxout-1,nyout-1)

          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var2d(i,j) = roufns(i+ia,j+ja)
            END DO
          END DO
          CALL netwrt2d(netid,packed,0,'ROUFNS',var2d,nxout-1,nyout-1)

          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var2d(i,j) = veg(i+ia,j+ja)
            END DO
          END DO
          CALL netwrt2d(netid,packed,0,'VEG',var2d,nxout-1,nyout-1)

        END IF
      END IF

      IF ( grdbas == 1 ) GOTO 3333

      CALL netwrtTime(netid,itime,'Time',curtim)

      IF (varout == 1) THEN

        IF (totout == 0) THEN

          DO k=1,nz-1
            DO j=1,nyout-1
              DO i=1,nxout
                var3du(i,j,k)=u(i+ia,j+ja,k)-ubar(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'UPRT',var3du,nxout,nyout-1,nz-1)

          DO k=1,nz-1
            DO j=1,nyout
              DO i=1,nxout-1
                var3dv(i,j,k)=v(i+ia,j+ja,k)-vbar(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'VPRT',var3dv,nxout-1,nyout,nz-1)

          DO k = 1,nz
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dw(i,j,k) = w(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'WPRT',var3dw,nxout-1,nyout-1,nz)

          DO k = 1,nz-1
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dt(i,j,k) = ptprt(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'PTPRT',var3dt,nxout-1,nyout-1,nz-1)

          DO k = 1,nz-1
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dt(i,j,k) = pprt(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'PPRT',var3dt,nxout-1,nyout-1,nz-1)

        ELSE

          DO k = 1,nz-1
            DO j = 1,nyout-1
              DO i = 1,nxout
                var3du(i,j,k) = u(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'U',var3du,nxout,nyout-1,nz-1)

          DO k = 1,nz-1
            DO j = 1,nyout
              DO i = 1,nxout-1
                var3dv(i,j,k) = v(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'V',var3dv,nxout-1,nyout,nz-1)

          DO k = 1,nz
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dw(i,j,k) = w(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'W',var3dw,nxout-1,nyout-1,nz)

          DO k=1,nz-1
            DO j=1,nyout-1
              DO i=1,nxout-1
                var3dt(i,j,k) = ptbar(i+ia,j+ja,k) + ptprt(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'PT',var3dt,nxout-1,nyout-1,nz-1)

          DO k=1,nz-1
            DO j=1,nyout-1
              DO i=1,nxout-1
                var3dt(i,j,k) = pbar(i+ia,j+ja,k) + pprt(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'P',var3dt,nxout-1,nyout-1,nz-1)

        END IF     ! totout
      END IF    ! varout

      IF (mstout == 1) THEN

        IF (totout == 0) THEN

          DO k=1,nz-1
            DO j=1,nyout-1
              DO i=1,nxout-1
                var3dt(i,j,k)=qv(i+ia,j+ja,k)-qvbar(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'QVPRT',var3dt,nxout-1,nyout-1,nz-1)

        ELSE

          DO k = 1,nz-1
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dt(i,j,k) = qv(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,'QV',var3dt,nxout-1,nyout-1,nz-1)

        END IF

        DO nq = 1, nscalar

          DO k = 1,nz-1
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dt(i,j,k) = qscalar(i+ia,j+ja,k,nq)
              END DO
            END DO
          END DO
          CALL netwrt3d(netid,packed,itime,upcase(qnames(nq)),var3dt,nxout-1,nyout-1,nz-1)
        END DO

        IF (rainout == 1) THEN

          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var2d(i,j) = raing(i+ia,j+ja)
            END DO
          END DO
          CALL netwrt2d(netid,packed,itime,'RAING',var2d,nxout-1,nyout-1)

          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var2d(i,j) = rainc(i+ia,j+ja)
            END DO
          END DO
          CALL netwrt2d(netid,packed,itime,'RAINC',var2d,nxout-1,nyout-1)

        END IF

        IF (prcout == 1) THEN

          DO k = 1,4
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3dw(i,j,k) = prcrate(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt2d(netid,packed,itime,'PRCRATE1',var3dw(:,:,1),nxout-1,nyout-1)
          CALL netwrt2d(netid,packed,itime,'PRCRATE2',var3dw(:,:,2),nxout-1,nyout-1)
          CALL netwrt2d(netid,packed,itime,'PRCRATE3',var3dw(:,:,3),nxout-1,nyout-1)
          CALL netwrt2d(netid,packed,itime,'PRCRATE4',var3dw(:,:,4),nxout-1,nyout-1)

        END IF

      END IF      ! mstout == 1

      IF (tkeout == 1) THEN

        DO k = 1,nz-1
          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var3dt(i,j,k) = tke(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'TKE',var3dt,nxout-1,nyout-1,nz-1)

      END IF

      IF (trbout == 1) THEN

        DO k = 1,nz-1
          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var3dt(i,j,k) = kmh(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'KMH',var3dt,nxout-1,nyout-1,nz-1)

        DO k = 1,nz-1
          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var3dt(i,j,k) = kmv(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'KMV',var3dt,nxout-1,nyout-1,nz-1)

      END IF  ! trbout

      IF (sfcout == 1) THEN

        DO is = 0,nstyps
          DO k = 1,nzsoil
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var4d(i,j,k,is+1) = tsoil(i+ia,j+ja,k,is)
              END DO
            END DO
          END DO
        END DO
        CALL netwrt4d(netid,packed,itime,'TSOIL',var4d,nxout-1,nyout-1,nzsoil,nstyps+1)

        DO is = 0,nstyps
          DO k = 1,nzsoil
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var4d(i,j,k,is+1) = qsoil(i+ia,j+ja,k,is)
              END DO
            END DO
          END DO
        END DO
        CALL netwrt4d(netid,packed,itime,'QSOIL',var4d,nxout-1,nyout-1,nzsoil,nstyps+1)

        DO is = 0,nstyps
          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var3dw(i,j,is+1) = wetcanp(i+ia,j+ja,is)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'WETCANP',var3dw,nxout-1,nyout-1,nstyps+1)

        IF (snowout == 1) THEN

          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var2d(i,j) = snowdpth(i+ia,j+ja)
            END DO
          END DO
          CALL netwrt2d(netid,packed,itime,'SNOWDPTH',var2d,nxout-1,nyout-1)

        END IF

      END IF  ! sfcout

      IF (radout == 1) THEN

        DO k = 1,nz-1
          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var3dt(i,j,k) = radfrc(i+ia,j+ja,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(netid,packed,itime,'RADFRC',var3dt,nxout-1,nyout-1,nz-1)

        DO j = 1,nyout-1
          DO i = 1,nxout-1
            var2d(i,j) = radsw(i+ia,j+ja)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'RADSW',var2d,nxout-1,nyout-1)

        DO j = 1,nyout-1
          DO i = 1,nxout-1
            var2d(i,j) = rnflx(i+ia,j+ja)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'RNFLX',var2d,nxout-1,nyout-1)

        DO j = 1,nyout-1
          DO i = 1,nxout-1
            var2d(i,j) = radswnet(i+ia,j+ja)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'RADSWNET',var2d,nxout-1,nyout-1)

        DO j = 1,nyout-1
          DO i = 1,nxout-1
            var2d(i,j) = radlwin(i+ia,j+ja)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'RADLWIN',var2d,nxout-1,nyout-1)

      END IF

      IF (flxout == 1) THEN

        DO j = 1,nyout-1
          DO i = 1,nxout-1
            var2d(i,j) = usflx(i+ia,j+ja)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'USFLX',var2d,nxout-1,nyout-1)

        DO j = 1,nyout-1
          DO i = 1,nxout-1
            var2d(i,j) = vsflx(i+ia,j+ja)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'VSFLX',var2d,nxout-1,nyout-1)

        DO j = 1,nyout-1
          DO i = 1,nxout-1
            var2d(i,j) = ptsflx(i+ia,j+ja)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'PTSFLX',var2d,nxout-1,nyout-1)

        DO j = 1,nyout-1
          DO i = 1,nxout-1
            var2d(i,j) = qvsflx(i+ia,j+ja)
          END DO
        END DO
        CALL netwrt2d(netid,packed,itime,'QVSFLX',var2d,nxout-1,nyout-1)

      END IF

      3333 CONTINUE
    END DO
  END DO

  DEALLOCATE(var2d,  var2di)
  DEALLOCATE(var3du, var3dv, var3dw, var3dt, var3di)
  DEALLOCATE(var4d)

  RETURN
END SUBROUTINE netsplitdump
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE netopen                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netopen(filename,fmode,nout)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!   Open a NetCDF file according to fmode
!
!     'C': Create a new NetCDF file
!     'W': Write to an exist NetCDF file
!     'R': Read from an exist NetCDF file
!
!------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: filename
  CHARACTER(LEN=1), INTENT(IN)  :: fmode
  INTEGER,          INTENT(OUT) :: nout

  INCLUDE    'netcdf.inc'

  INTEGER :: istatus
  LOGICAL :: fexists
  LOGICAL :: LargeFile

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  LargeFile = .TRUE.

  SELECT CASE (fmode)

    CASE ('C','W')

      IF (LargeFile) THEN
        istatus = NF_CREATE(TRIM(filename),IOR(NF_CLOBBER,NF_64BIT_OFFSET),&
                        nout)                                    ! CDF2
      ELSE
        istatus = NF_CREATE(TRIM(filename),NF_CLOBBER,nout)      ! CDF1
      END IF
      CALL net_check_error(istatus,'netopen')

    CASE ('R')

      INQUIRE(FILE = TRIM(filename), EXIST = fexists)
      IF (fexists) THEN
        istatus = NF_OPEN(TRIM(filename),NF_NOWRITE,nout)
        CALL net_check_error(istatus,'netopen')
      ELSE
        WRITE(6,'(2a)') 'File not found: ', filename
        istatus = -1
      END IF

    CASE ('A')  ! file should be in define mode after this call

      INQUIRE(FILE = TRIM(filename), EXIST = fexists)
      IF (fexists) THEN
        istatus = NF_OPEN(TRIM(filename),NF_WRITE,nout)
        CALL net_check_error(istatus,'netopen')

        istatus = NF_REDEF(nout)
        CALL net_check_error(istatus,'netopen')
      ELSE
        WRITE(6,'(2a)') 'File not found: ', filename
        istatus = -1
      END IF

    CASE DEFAULT

      WRITE(6,*) 'Wrong file mode: ',fmode,'.'
      istatus = -1

  END SELECT

  IF (istatus < 0) nout = istatus !CALL arpsstop('**** ERROR in netopen ****',1)

  RETURN
END SUBROUTINE netopen
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE net_check_error           ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_check_error(ierr,sub_name)

  IMPLICIT NONE

  INTEGER,          INTENT(IN) :: ierr
  CHARACTER(LEN=*), INTENT(IN) :: sub_name

  CHARACTER(LEN=80) :: errmsg

  INCLUDE 'netcdf.inc'

  IF(ierr /= NF_NOERR) THEN
    errmsg = NF_STRERROR(ierr)
    WRITE(6,'(/2a)') 'NetCDF error: ',errmsg
    WRITE(6,'(3a/)' ) 'Program stopped while calling "', sub_name,'".'
    STOP
  END IF

  RETURN
END SUBROUTINE net_check_error
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE netclose                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netclose(nchout)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!   Close the NetCDF file.
!
!------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: nchout

  INCLUDE 'netcdf.inc'

  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  istatus = NF_CLOSE(nchout)
  CALL net_check_error(istatus,'netclose')

  RETURN
END SUBROUTINE netclose
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_define_dimension           ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_define_dimension(ncid,grdbas,nx,ny,nz,nzsoil,nstyps)
!
!-----------------------------------------------------------------------
!
!  Define dimensions for ARPS History file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(IN)  :: grdbas
  INTEGER, INTENT(IN)  :: nx
  INTEGER, INTENT(IN)  :: ny
  INTEGER, INTENT(IN)  :: nz
  INTEGER, INTENT(IN)  :: nzsoil
  INTEGER, INTENT(IN)  :: nstyps

  INCLUDE    'netcdf.inc'

  INTEGER  :: istatus
  INTEGER  :: dimunlim_id
  INTEGER  :: dimwe_id, dimwes_id, dimsn_id, dimsns_id
  INTEGER  :: dimbt_id, dimbts_id
  INTEGER  :: dimsoil_id,  dimn_id,  dimns_id

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ! define dimensions
  IF (grdbas /= 1) THEN
    istatus = NF_DEF_DIM(ncid,'Time',NF_UNLIMITED,dimunlim_id)
    CALL net_check_error(istatus,'net_define_dimension')
  END IF

  istatus = NF_DEF_DIM(ncid,'x',nx-1,dimwe_id)
  CALL net_check_error(istatus,'net_define_dimension')

  istatus = NF_DEF_DIM(ncid,'y',ny-1,dimsn_id)
  CALL net_check_error(istatus,'net_define_dimension')

  istatus = NF_DEF_DIM(ncid,'z',nz-1,dimbt_id)
  CALL net_check_error(istatus,'net_define_dimension')

  istatus = NF_DEF_DIM(ncid,'x_stag',nx,dimwes_id)
  CALL net_check_error(istatus,'net_define_dimension')

  istatus = NF_DEF_DIM(ncid,'y_stag',ny,dimsns_id)
  CALL net_check_error(istatus,'net_define_dimension')

  istatus = NF_DEF_DIM(ncid,'z_stag',nz,dimbts_id)
  CALL net_check_error(istatus,'net_define_dimension')

  istatus = NF_DEF_DIM(ncid,'zsoil',nzsoil,dimsoil_id)
  CALL net_check_error(istatus,'net_define_dimension')

  istatus = NF_DEF_DIM(ncid,'nstyp',nstyps,dimn_id)
  CALL net_check_error(istatus,'net_define_dimension')

  IF (grdbas /= 1) THEN
    istatus = NF_DEF_DIM(ncid,'nstyp_total',nstyps+1,dimns_id)
    CALL net_check_error(istatus,'net_define_dimension')
  END IF

  RETURN
END SUBROUTINE net_define_dimension
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_define_variables           ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_define_variables(ncid,packed,grdbas,tmunit,nx,ny,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!     Define ARPS history file attributes and variables. After this call
!     The netCDF file should be in DATA mode.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (08/10/2004)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  USE arps_netio_metadata

  IMPLICIT NONE
  INTEGER, INTENT(IN)          :: ncid
  INTEGER, INTENT(IN)          :: packed         ! may support pack latter
  INTEGER, INTENT(IN)          :: grdbas
  INTEGER, INTENT(IN)          :: nx,ny          ! must be global domain size
  CHARACTER(LEN=*), INTENT(IN) :: tmunit
  INTEGER, INTENT(OUT)         :: istatus

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'           ! Grid & map parameters.
  INCLUDE 'phycst.inc'

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER           :: lenstr
  INTEGER           :: i,j,k,n
  INTEGER           :: nq
  INTEGER           :: varid
  INTEGER           :: dimns_id,dimunlim_id
  INTEGER           :: dimx_id, dimy_id, dimz_id, dimsoil_id,dimn_id
  INTEGER           :: dimxs_id,dimys_id,dimzs_id

  INTEGER           :: oldfillmode

  CHARACTER(LEN=80) :: tmpstr,tstr

  CHARACTER(LEN=4) :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (grdbas == 1) THEN
    tmpstr   = 'ARPS 5.3 grid & base (time independent) data'
  ELSE
    tmpstr   = 'ARPS 5.3 history dump (time dependent)'
  END IF

  CALL netwrt_general_att(ncid,packed,tmpstr,nx,ny,dx,dy,mapproj,sclfct,&
                          trulat1,trulat2,trulon,ctrlat,ctrlon,istatus)

!-----------------------------------------------------------------------
!
! Define specific global attributes for ARPS history files
!
!-----------------------------------------------------------------------

  lenstr  = LEN_TRIM(runname)
  istatus = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'RUNNAME',lenstr,runname)
  istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'nocmnt',NF_INT,1,nocmnt)
  tmpstr(:)  = ' '
  DO n = 1, nocmnt
    WRITE(tmpstr,'(a,I2.2)') 'cmnt',n
    lenstr  = LEN_TRIM(cmnt(n))
    istatus = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,TRIM(tmpstr),lenstr,cmnt(n))
  END DO

  !
  ! Date & time
  !
  tmpstr(:) = ' '
  WRITE(tmpstr,'(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2,a,I2.2)')             &
                  year,'-', month,'-',day,'_',hour,':',minute,':',second
  lenstr  = LEN_TRIM(tmpstr)
  istatus = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'INITIAL_TIME',lenstr,tmpstr)

  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'TSTOP',  NF_FLOAT,1,tstop)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'THISDMP',NF_FLOAT,1,thisdmp)

  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'LATITUD',   NF_FLOAT,1,latitud)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'XGRDORG',   NF_FLOAT,1,xgrdorg)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'YGRDORG',   NF_FLOAT,1,ygrdorg)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'UMOVE',     NF_FLOAT,1,umove)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'VMOVE',     NF_FLOAT,1,vmove)

  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'NTCLOUD',    NF_FLOAT,1,ntcloud)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'N0RAIN',    NF_FLOAT,1,n0rain)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'N0SNOW',    NF_FLOAT,1,n0snow)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'N0GRPL',    NF_FLOAT,1,n0grpl)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'N0HAIL',    NF_FLOAT,1,n0hail)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'RHOICE',   NF_FLOAT,1,rhoice)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'RHOSNOW',   NF_FLOAT,1,rhosnow)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'RHOGRPL',   NF_FLOAT,1,rhogrpl)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'RHOHAIL',   NF_FLOAT,1,rhohail)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'ALPHARAIN',   NF_FLOAT,1,alpharain)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'ALPHAICE',   NF_FLOAT,1,alphaice)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'ALPHASNOW',   NF_FLOAT,1,alphasnow)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'ALPHAGRPL',   NF_FLOAT,1,alphagrpl)
  istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'ALPHAHAIL',   NF_FLOAT,1,alphahail)

  !
  ! Flags
  !
  IF( grdbas == 1 ) THEN
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'GRDFLG', NF_INT,1,1)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'BASFLG', NF_INT,1,1)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'VARFLG', NF_INT,1,0)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'MSTFLG', NF_INT,1,mstout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'ICEFLG', NF_INT,1,0)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'TRBFLG', NF_INT,1,0)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'SFCFLG', NF_INT,1,0)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'RAINFLG',NF_INT,1,0)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'LANDFLG',NF_INT,1,landout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'TOTFLG', NF_INT,1,totout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'TKEFLG', NF_INT,1,0)
  ELSE
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'GRDFLG', NF_INT,1,grdout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'BASFLG', NF_INT,1,basout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'VARFLG', NF_INT,1,varout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'MSTFLG', NF_INT,1,mstout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'ICEFLG', NF_INT,1,iceout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'TRBFLG', NF_INT,1,trbout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'SFCFLG', NF_INT,1,sfcout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'RAINFLG',NF_INT,1,rainout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'LANDFLG',NF_INT,1,landout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'TOTFLG', NF_INT,1,totout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'TKEFLG', NF_INT,1,tkeout)

    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'nscalar',NF_INT,1,nscalar)

    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_QC',NF_INT,1,P_QC)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_QR',NF_INT,1,P_QR)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_QI',NF_INT,1,P_QI)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_QS',NF_INT,1,P_QS)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_QG',NF_INT,1,P_QG)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_QH',NF_INT,1,P_QH)

    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_NC',NF_INT,1,P_NC)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_NR',NF_INT,1,P_NR)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_NI',NF_INT,1,P_NI)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_NS',NF_INT,1,P_NS)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_NG',NF_INT,1,P_NG)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_NH',NF_INT,1,P_NH)

    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_ZR',NF_INT,1,P_ZR)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_ZI',NF_INT,1,P_ZI)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_ZS',NF_INT,1,P_ZS)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_ZG',NF_INT,1,P_ZG)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_ZH',NF_INT,1,P_ZH)

    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'P_CC',NF_INT,1,P_CC)
  END IF

  IF ( totout == 1 ) THEN
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'PRCFLG', NF_INT,1,prcout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'RADFLG', NF_INT,1,radout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'FLXFLG', NF_INT,1,flxout)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'SNOWFLG',NF_INT,1,snowout)
  END IF

  ! do not fill, will set values explicitly later. Improve performance

  istatus = NF_SET_FILL(ncid,NF_NOFILL,oldfillmode)

!-----------------------------------------------------------------------
!
! Define variable arrays
!
!-----------------------------------------------------------------------

  !
  ! Get dimension IDs
  !
  istatus = NF_INQ_DIMID(ncid,'x_stag', dimx_id)
  istatus = NF_INQ_DIMID(ncid,'y_stag', dimy_id)
  istatus = NF_INQ_DIMID(ncid,'z_stag', dimz_id)
  istatus = NF_INQ_DIMID(ncid,'x',      dimxs_id)
  istatus = NF_INQ_DIMID(ncid,'y',      dimys_id)
  istatus = NF_INQ_DIMID(ncid,'z',      dimzs_id)
  istatus = NF_INQ_DIMID(ncid,'zsoil',  dimsoil_id)
  istatus = NF_INQ_DIMID(ncid,'nstyp',  dimn_id)

  IF(grdbas /= 1)  THEN

    istatus = NF_INQ_DIMID(ncid, 'Time',            dimunlim_id)
    CALL net_check_error(istatus,'net_define_variabls:Time')
    istatus = NF_INQ_DIMID(ncid, 'nstyp_total',dimns_id)
    CALL net_check_error(istatus,'net_define_variabls:nstyp_total')

  END IF

  IF(grdout == 1 .OR. grdbas == 1 ) THEN
    !
    ! Grid variables
    !
    istatus = NF_DEF_VAR(ncid,'x_stag',NF_FLOAT,1,(/dimx_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%x_stag)

    istatus = NF_DEF_VAR(ncid,'y_stag',NF_FLOAT,1,(/dimy_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%y_stag)

    istatus = NF_DEF_VAR(ncid,'z_stag',NF_FLOAT,1,(/dimz_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%z_stag)

    istatus = NF_DEF_VAR(ncid,'ZP',NF_FLOAT,3,(/dimxs_id,dimys_id,dimz_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%zp)

    istatus = NF_DEF_VAR(ncid,'ZPSOIL',NF_FLOAT,3,(/dimxs_id,dimys_id,dimsoil_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%zpsoil)
  END IF

  IF(basout == 1 .OR. grdbas == 1 ) THEN
    !
    ! Base state variables
    !
    istatus = NF_DEF_VAR(ncid,'UBAR',NF_FLOAT,3,(/dimx_id,dimys_id,dimzs_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%ubar)

    istatus = NF_DEF_VAR(ncid,'VBAR',NF_FLOAT,3,(/dimxs_id,dimy_id,dimzs_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%vbar)

    istatus = NF_DEF_VAR(ncid,'WBAR',NF_FLOAT,3,(/dimxs_id,dimys_id,dimz_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%wbar)

    istatus = NF_DEF_VAR(ncid,'PTBAR',NF_FLOAT,3,(/dimxs_id,dimys_id,dimzs_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%ptbar)

    istatus = NF_DEF_VAR(ncid,'PBAR',NF_FLOAT,3,(/dimxs_id,dimys_id,dimzs_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%pbar)

    IF(mstout == 1) THEN
      istatus = NF_DEF_VAR(ncid,'QVBAR',NF_FLOAT,3,(/dimxs_id,dimys_id,dimzs_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%qvbar)
    END IF

    IF(landout == 1) THEN

      istatus = NF_DEF_VAR(ncid,'SOILTYP',NF_INT,3,(/dimxs_id,dimys_id,dimn_id/),varid)
      CALL net_define_var_meta(ncid,varid,'INT',arpsmeta%soiltyp)

      istatus = NF_DEF_VAR(ncid,'STYPFRCT',NF_FLOAT,3,(/dimxs_id,dimys_id,dimn_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%stypfrct)

      istatus = NF_DEF_VAR(ncid,'VEGTYP',NF_INT,2,(/dimxs_id,dimys_id/),varid)
      CALL net_define_var_meta(ncid,varid,'INT',arpsmeta%vegtyp)

      istatus = NF_DEF_VAR(ncid,'LAI',NF_FLOAT,2,(/dimxs_id,dimys_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%lai)

      istatus = NF_DEF_VAR(ncid,'ROUFNS',NF_FLOAT,2,(/dimxs_id,dimys_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%roufns)

      istatus = NF_DEF_VAR(ncid,'VEG',NF_FLOAT,2,(/dimxs_id,dimys_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%veg)

    END IF
  END IF

  IF ( grdbas == 1 ) GOTO 444    ! For grid & base file, this is all of the data to be written

  tmpstr(:) = ' '
  WRITE(tmpstr,'(I4.4,a,I2.2,a,I2.2,1x,I2.2,a,I2.2,a,I2.2)')            &
                  year,'-', month,'-',day,hour,':',minute,':',second
  WRITE(tstr,'(3a)') TRIM(tmunit), ' since ',TRIM(tmpstr)

  istatus = NF_DEF_VAR(ncid,'Time',NF_FLOAT,1,(/dimunlim_id/),varid)

  arpsmeta%Time%units = tstr
  CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%Time)

  IF(varout == 1) THEN

    IF ( totout == 0 ) THEN

      istatus = NF_DEF_VAR(ncid,'UPRT',NF_FLOAT,4,(/dimx_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%uprt)

      istatus = NF_DEF_VAR(ncid,'VPRT',NF_FLOAT,4,(/dimxs_id,dimy_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%vprt)

      istatus = NF_DEF_VAR(ncid,'WPRT',NF_FLOAT,4,(/dimxs_id,dimys_id,dimz_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%wprt)

      istatus = NF_DEF_VAR(ncid,'PTPRT',NF_FLOAT,4,(/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%ptprt)

      istatus = NF_DEF_VAR(ncid,'PPRT',NF_FLOAT,4,(/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%pprt)

    ELSE

      istatus = NF_DEF_VAR(ncid,'U',NF_FLOAT,4,(/dimx_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%u)

      istatus = NF_DEF_VAR(ncid,'V',NF_FLOAT,4,(/dimxs_id,dimy_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%v)

      istatus = NF_DEF_VAR(ncid,'W',NF_FLOAT,4,(/dimxs_id,dimys_id,dimz_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%w)

      istatus = NF_DEF_VAR(ncid,'PT',NF_FLOAT,4,(/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%pt)

      istatus = NF_DEF_VAR(ncid,'P',NF_FLOAT,4,(/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%p)

    END IF
  END IF

  IF(mstout == 1) THEN

    IF( totout == 0 ) THEN
      istatus = NF_DEF_VAR(ncid,'QVPRT',NF_FLOAT,4,(/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%qvprt)
    ELSE
      istatus = NF_DEF_VAR(ncid,'QV',NF_FLOAT,4,(/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%qv)
    END IF

    DO nq = 1,nscalar
      istatus = NF_DEF_VAR(ncid,upcase(qnames(nq)),NF_FLOAT,4,    &
                       (/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%qscalar(nq))
    END DO

    IF(rainout == 1) THEN

      istatus = NF_DEF_VAR(ncid,'RAING',NF_FLOAT,3,(/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%raing)

      istatus = NF_DEF_VAR(ncid,'RAINC',NF_FLOAT,3,(/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%rainc)

    END IF

    IF ( prcout == 1 ) THEN

      istatus = NF_DEF_VAR(ncid,'PRCRATE1',NF_FLOAT,3,(/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%prcrate1)

      istatus = NF_DEF_VAR(ncid,'PRCRATE2',NF_FLOAT,3,(/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%prcrate2)

      istatus = NF_DEF_VAR(ncid,'PRCRATE3',NF_FLOAT,3,(/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%prcrate3)

      istatus = NF_DEF_VAR(ncid,'PRCRATE4',NF_FLOAT,3,(/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%prcrate4)

    END IF

  END IF       ! mstout == 1

  IF( tkeout == 1 ) THEN

    istatus = NF_DEF_VAR(ncid,'TKE',NF_FLOAT,4,(/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%tke)

  END IF

  IF( trbout == 1 ) THEN

    istatus = NF_DEF_VAR(ncid,'KMH',NF_FLOAT,4,(/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%kmh)

    istatus = NF_DEF_VAR(ncid,'KMV',NF_FLOAT,4,(/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%kmv)

  END IF

  IF( sfcout == 1) THEN

    istatus = NF_DEF_VAR(ncid,'TSOIL',NF_FLOAT,5,                       &
            (/dimxs_id,dimys_id,dimsoil_id,dimns_id,dimunlim_id/),varid)
    CALL net_check_error(istatus,'net_define_variabls:TSOIL')
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%tsoil)

    istatus = NF_DEF_VAR(ncid,'QSOIL',NF_FLOAT,5,                       &
            (/dimxs_id,dimys_id,dimsoil_id,dimns_id,dimunlim_id/),varid)
    CALL net_check_error(istatus,'net_define_variabls:QSOIL')
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%qsoil)

    istatus = NF_DEF_VAR(ncid,'WETCANP',NF_FLOAT,4,                     &
                       (/dimxs_id,dimys_id,dimns_id,dimunlim_id/),varid)
    CALL net_check_error(istatus,'net_define_variabls:WETCANP')
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%wetcanp)

    IF (snowout == 1) THEN
      istatus = NF_DEF_VAR(ncid,'SNOWDPTH',NF_FLOAT,3,                  &
                                (/dimxs_id,dimys_id,dimunlim_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%snowdpth)
    END IF

  END IF           ! sfcout == 1

  IF( radout == 1 ) THEN

    istatus = NF_DEF_VAR(ncid,'RADFRC',NF_FLOAT,4,                      &
                         (/dimxs_id,dimys_id,dimzs_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%radfrc)

    istatus = NF_DEF_VAR(ncid,'RADSW',NF_FLOAT,3,                       &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%radsw)

    istatus = NF_DEF_VAR(ncid,'RNFLX',NF_FLOAT,3,                       &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%rnflx)

    istatus = NF_DEF_VAR(ncid,'RADSWNET',NF_FLOAT,3,                    &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%radswnet)

    istatus = NF_DEF_VAR(ncid,'RADLWIN',NF_FLOAT,3,                     &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%radlwin)

  END IF

  IF( flxout == 1 ) THEN

    istatus = NF_DEF_VAR(ncid,'USFLX',NF_FLOAT,3,                       &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%usflx)

    istatus = NF_DEF_VAR(ncid,'VSFLX',NF_FLOAT,3,                       &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%vsflx)

    istatus = NF_DEF_VAR(ncid,'PTSFLX',NF_FLOAT,3,                      &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%ptsflx)

    istatus = NF_DEF_VAR(ncid,'QVSFLX',NF_FLOAT,3,                       &
                                 (/dimxs_id,dimys_id,dimunlim_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',arpsmeta%qvsflx)

  END IF

!-----------------------------------------------------------------------
!
! End NetCDF file DEFINE mode
!
!-----------------------------------------------------------------------

  444 CONTINUE

  istatus = NF_ENDDEF(ncid)
  CALL net_check_error(istatus,'net_define_variabls')

  RETURN
END SUBROUTINE net_define_variables
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE netwrtTime               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netwrtTime(nout,itime,varname,var1d)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 1D vector to the output file.
!
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)          :: nout   ! output channel,  NetCDF id
  INTEGER, INTENT(IN)          :: itime
  REAL,    INTENT(IN)          :: var1d
  CHARACTER(LEN=*), INTENT(IN) :: varname

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid, istatus

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing data valid time ', varname

  !
  ! get variable id and dimension length
  !
  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in netwrtTime.')

  !
  ! Write data
  !
  istatus = NF_PUT_VARA_REAL(nout,varid,(/itime/),(/1/),var1d)
  CALL net_check_error(istatus,'NF_PUT_VARA_REAL in netwrtTime')

  WRITE(6,'(14x,a)') '=== DONE ==='

  RETURN
END SUBROUTINE netwrtTime
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE netwrt1d                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netwrt1d(nout,packed,itime,varname,var1d,ndim)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 1D vector to the output file.
!
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)          :: nout   ! output channel,  NetCDF id
  INTEGER, INTENT(IN)          :: packed
  INTEGER, INTENT(IN)          :: itime
  INTEGER, INTENT(IN)          :: ndim
  REAL,    INTENT(IN)          :: var1d(ndim)
  CHARACTER(LEN=*), INTENT(IN) :: varname

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid, istatus
  INTEGER :: dim_ids(2)
  INTEGER :: dimlens(2)

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 1D NetCDF variable ', varname

  !
  ! get variable id and dimension length
  !
  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in write1d.')
  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)
  istatus = NF_INQ_DIMLEN(nout,dim_ids(1),dimlens(1))

  !
  ! check dimension
  !
  IF(dimlens(1) /= ndim) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndim,               &
                   ', Defined dimension in file = ',dimlens(1)
    STOP
  END IF

  !
  ! Write data
  !
  IF (itime > 0) THEN
    istatus = NF_PUT_VARA_REAL(nout,varid,(/1,itime/),(/ndim,1/),var1d)
  ELSE
    istatus = NF_PUT_VARA_REAL(nout,varid,(/1/),(/ndim/),var1d)
  END IF
  CALL net_check_error(istatus,'NF_PUT_VARA_REAL in netwrt1d')

  WRITE(6,'(a)') '         === DONE ==='

  RETURN
END SUBROUTINE netwrt1d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE netwrt2d                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netwrt2d(nout,packed,itime,varname,var2d,ndimx,ndimy)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 2D array to the output file.
!
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nout
  INTEGER, INTENT(IN) :: packed
  INTEGER, INTENT(IN) :: itime             ! = 0  no unlimited dimension
                                           ! > 0  record No.
  INTEGER, INTENT(IN) :: ndimx,ndimy
  REAL,    INTENT(IN) :: var2d(ndimx,ndimy)
  CHARACTER(LEN=*), INTENT(IN) :: varname

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid, istatus
  INTEGER :: dim_ids(3)
  INTEGER :: dimlens(3)

  CHARACTER(LEN=20) :: fmtstr

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 2D NetCDF variable ', varname

! get variable id and dimension length

  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in write2d.')

  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)

  istatus = NF_INQ_DIMLEN(nout,dim_ids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(2),dimlens(2))
  IF(itime > 0)         &
  istatus = NF_INQ_DIMLEN(nout,dim_ids(3),dimlens(3))   ! unlimit dimension

  IF(dimlens(1) /= ndimx) THEN
    WRITE(6,'(/a)') ' ERROR: Mismatched dimension size in X direction.'
    WRITE(6,*) 'Input X dimension = ',ndimx, ' Data X dimension =',dimlens(1)
    STOP
  END IF

  IF(dimlens(2) /= ndimy) THEN
    WRITE(6,'(/a)') ' ERROR:Mismatched dimension size in Y direction.'
    WRITE(6,*) 'Input Y dimension = ',ndimx, ' Data Y dimension =',dimlens(2)
    STOP
  END IF

! Write data

  IF (itime > 0) THEN
    istatus = NF_PUT_VARA_REAL(nout,varid,(/1,1,itime/),                &
                                          (/ndimx,ndimy,1/),var2d)
  ELSE
    istatus = NF_PUT_VARA_REAL(nout,varid,(/1,1/),(/ndimx,ndimy/),var2d)
  END IF
  CALL net_check_error(istatus,'NF_PUT_VARA_REAL in netwrt2d')

  WRITE(fmtstr,'(a,I2.2,a)') '(',15-LEN_TRIM(varname),'x,a)'
  WRITE(6,FMT=fmtstr) '=== DONE ==='

  RETURN
END SUBROUTINE netwrt2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE netwrt2di                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netwrt2di(nout,packed,itime,varname,var2d,ndimx,ndimy)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 2D array to the output file.
!
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nout
  INTEGER, INTENT(IN) :: packed
  INTEGER, INTENT(IN) :: itime
  INTEGER, INTENT(IN) :: ndimx,ndimy
  INTEGER, INTENT(IN) :: var2d(ndimx,ndimy)

  CHARACTER(LEN=*), INTENT(IN) :: varname

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid, istatus
  INTEGER :: dim_ids(3)
  INTEGER :: dimlens(3)

  CHARACTER(LEN=20) :: fmtstr

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,FMT='(2a)',ADVANCE='NO')     &
                        '  Writing 2D integer variable ', varname

  !
  ! get variable id and dimension length
  !
  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in write2di.')
  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)
  istatus = NF_INQ_DIMLEN(nout,dim_ids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(2),dimlens(2))
  IF (itime > 0)                      &              ! unlimit dimension
  istatus = NF_INQ_DIMLEN(nout,dim_ids(3),dimlens(3))

  !
  ! Some checks to confirm it is the right variable and dimensions
  !
  IF(dimlens(1) /= ndimx) THEN
    WRITE(6,'(/a)') 'Mismatched dimension size in X direction.'
    WRITE(6,*) ' Input X dimension = ',ndimx,                           &
               ' Data X dimension = ', dimlens(1)
    STOP
  END IF

  IF(dimlens(2) /= ndimy) THEN
    WRITE(6,'(/a)') 'Mismatched dimension size in Y direction.'
    WRITE(6,*) ' Input Y dimension = ',ndimy,                           &
               ' Data Y dimension = ',dimlens(2)
    STOP
  END IF

  !
  ! Write data
  !
  IF (itime > 0) THEN
    istatus = NF_PUT_VARA_INT(nout,varid,(/1,1,itime/),                 &
                                         (/ndimx,ndimy,1/),var2d)
  ELSE
    istatus = NF_PUT_VARA_INT(nout,varid,(/1,1/),(/ndimx,ndimy/),var2d)
  END IF
  CALL net_check_error(istatus,'NF_PUT_VARA_INT in netwrt2di.')

  WRITE(fmtstr,'(a,I2.2,a)') '(',14-LEN_TRIM(varname),'x,a)'
  WRITE(6,FMT=fmtstr) '=== DONE ==='

  RETURN
END SUBROUTINE netwrt2di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE netwrt3d                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netwrt3d(nout,packed,itime,varname,var3d,ndimx,ndimy,ndimz)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 3D array to the output file.
!
!------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nout
  INTEGER, INTENT(IN) :: packed
  INTEGER, INTENT(IN) :: itime
  INTEGER, INTENT(IN) :: ndimx,ndimy,ndimz
  REAL,    INTENT(IN) :: var3d(ndimx,ndimy,ndimz)

  CHARACTER(LEN=*), INTENT(IN) :: varname

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid, istatus
  INTEGER :: dim_ids(4)
  INTEGER :: dimlens(4)

  CHARACTER(LEN=20) :: fmtstr

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,FMT='(1x,2a)',ADVANCE='NO') ' Writing 3D NetCDF variable ', varname

  !
  ! get variable id
  !
  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in write3d.')

  !
  ! get dimension lengths and do some checks
  !
  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)
  istatus = NF_INQ_DIMLEN  (nout,dim_ids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN  (nout,dim_ids(2),dimlens(2))
  istatus = NF_INQ_DIMLEN  (nout,dim_ids(3),dimlens(3))
  IF(itime > 0)     &                                 ! unlimit dimension
  istatus = NF_INQ_DIMLEN(nout,dim_ids(4),dimlens(4))

  IF(dimlens(1) /= ndimx) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size in X direction.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndimx,              &
                   ', Defined dimension in file = ',dimlens(1)
    STOP
  END IF
  IF(dimlens(2) /= ndimy) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size in Y direction.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndimy,              &
                   ', Defined dimension in file = ',dimlens(2)
    STOP
  END IF
  IF(dimlens(3) /= ndimz) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size in the 3rd dimension.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndimz,              &
                   ', Defined dimension in file = ',dimlens(3)
    STOP
  END IF

  !
  ! Write data
  !
  IF (itime > 0 ) THEN             ! Actually, 4D array defined in file
    istatus = NF_PUT_VARA_REAL(nout,varid,(/1,1,1,itime/),              &
                              (/ndimx,ndimy,ndimz,1/),var3d)
  ELSE                             ! Just 3D array, no unlimited dim.
    istatus = NF_PUT_VARA_REAL(nout,varid,(/1,1,1/),                    &
                              (/ndimx,ndimy,ndimz/),var3d)
  END IF
  CALL net_check_error(istatus,'NF_PUT_VARA_REAL in netwrt3d')

  WRITE(fmtstr,'(a,I2.2,a)') '(',15-LEN_TRIM(varname),'x,a)'
  WRITE(6,FMT=fmtstr) '=== DONE ==='


  RETURN
END SUBROUTINE netwrt3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE netwrt3di                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netwrt3di(nout,packed,itime,varname,var3d,ndimx,ndimy,ndimz)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 3D integer array to the output file.
!
!------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nout
  INTEGER, INTENT(IN) :: packed
  INTEGER, INTENT(IN) :: itime
  INTEGER, INTENT(IN) :: ndimx,ndimy,ndimz
  INTEGER, INTENT(IN) :: var3d(ndimx,ndimy,ndimz)

  CHARACTER(LEN=*), INTENT(IN) :: varname

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid, istatus
  INTEGER :: dim_ids(4)
  INTEGER :: dimlens(4)

  CHARACTER(LEN=20) :: fmtstr

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,FMT='(1x,2a)',ADVANCE='NO')                                   &
                         ' Writing 3D integer variable ', varname
  !
  ! get variable id
  !
  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in netwrt3di.')

  !
  ! get dimension lengths and do some checks
  !
  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)
  istatus = NF_INQ_DIMLEN  (nout,dim_ids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN  (nout,dim_ids(2),dimlens(2))
  istatus = NF_INQ_DIMLEN  (nout,dim_ids(3),dimlens(3))
  IF(itime > 0)     &                                 ! unlimit dimension
  istatus = NF_INQ_DIMLEN(nout,dim_ids(4),dimlens(4))

  IF(dimlens(1) /= ndimx) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size in X direction.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndimx,              &
                   ', Defined dimension in file = ',dimlens(1)
    STOP
  END IF
  IF(dimlens(2) /= ndimy) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size in Y direction.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndimy,              &
                   ', Defined dimension in file = ',dimlens(2)
    STOP
  END IF
  IF(dimlens(3) /= ndimz) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size in the 3rd dimension.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndimz,              &
                   ', Defined dimension in file = ',dimlens(3)
    STOP
  END IF

  !
  ! Write data
  !
  IF (itime > 0 ) THEN             ! Actually, 4D array defined in file
    istatus = NF_PUT_VARA_INT(nout,varid,(/1,1,1,itime/),               &
                              (/ndimx,ndimy,ndimz,1/),var3d)
  ELSE                             ! Just 3D array, no unlimited dim.
    istatus = NF_PUT_VARA_INT(nout,varid,(/1,1,1/),                     &
                              (/ndimx,ndimy,ndimz/),var3d)
  END IF
  CALL net_check_error(istatus,'NF_PUT_VARA_INT in netwrt3di')

  WRITE(fmtstr,'(a,I2.2,a)') '(',14-LEN_TRIM(varname),'x,a)'
  WRITE(6,FMT=fmtstr) '=== DONE ==='

  RETURN
END SUBROUTINE netwrt3di
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE netwrt4d                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netwrt4d(nout,packed,itime,varname,var4d,                    &
                    ndimx,ndimy,ndimz,ndims)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 4D array to the output file.
!
!------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nout
  INTEGER, INTENT(IN) :: packed
  INTEGER, INTENT(IN) :: itime
  INTEGER, INTENT(IN) :: ndimx,ndimy,ndimz,ndims
  REAL,    INTENT(IN) :: var4d(ndimx,ndimy,ndimz,ndims)

  CHARACTER(LEN=*), INTENT(IN) :: varname

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid, istatus
  INTEGER :: dim_ids(5)
  INTEGER :: dimlens(5)

  CHARACTER(LEN=20) :: fmtstr

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,FMT='(1x,2a)',ADVANCE='NO') ' Writing 4D NetCDF variable ', varname

  !
  ! get variable id
  !
  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in netwrt4d.')

  !
  ! get dimension lengths and do some checks
  !
  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)
  istatus = NF_INQ_DIMLEN  (nout,dim_ids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN  (nout,dim_ids(2),dimlens(2))
  istatus = NF_INQ_DIMLEN  (nout,dim_ids(3),dimlens(3))
  istatus = NF_INQ_DIMLEN  (nout,dim_ids(4),dimlens(4))
  IF(itime > 0) istatus = NF_INQ_DIMLEN(nout,dim_ids(5),dimlens(5))
                                                     ! unlimit dimension

  IF(dimlens(1) /= ndimx) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size in X direction.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndimx,              &
                   ', Defined dimension in file = ',dimlens(1)
    STOP
  END IF
  IF(dimlens(2) /= ndimy) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size in Y direction.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndimy,              &
                   ', Defined dimension in file = ',dimlens(2)
    STOP
  END IF
  IF(dimlens(3) /= ndimz) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size in the 3rd dimension.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndimz,              &
                   ', Defined dimension in file = ',dimlens(3)
    STOP
  END IF
  IF(dimlens(4) /= ndims) THEN
    WRITE(6,'(1x,/a)') 'Mismatched dimension size in the 4th dimension.'
    WRITE(6,'(1x,2(a,I4))')   ' Input dimension = ',ndims,              &
                   ', Defined dimension in file = ',dimlens(4)
    STOP
  END IF

  !
  ! Write data
  !
  IF (itime > 0 ) THEN             ! Actually, 5D array defined in file
    istatus = NF_PUT_VARA_REAL(nout,varid,(/1,1,1,1,itime/),            &
                              (/ndimx,ndimy,ndimz,ndims,1/),var4d)
  ELSE                             ! Just 4D array, no unlimited dim.
    istatus = NF_PUT_VARA_REAL(nout,varid,(/1,1,1,1/),                  &
                              (/ndimx,ndimy,ndimz,ndims/),var4d)
  END IF
  CALL net_check_error(istatus,'NF_PUT_VARA_REAL in netwrt4d')

  WRITE(fmtstr,'(a,I2.2,a)') '(',15-LEN_TRIM(varname),'x,a)'
  WRITE(6,FMT=fmtstr) '=== DONE ==='

  RETURN
END SUBROUTINE netwrt4d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE netreaddims                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netreaddims(ncid,nxout,nyout,nzout,nzsoilout,          &
                       nstypsout,istatus)

!------------------------------------------------------------------------
!
! PURPOSE:
!
!   Read dimension parameters from NetCDF output file.
!
!------------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid

  INTEGER, INTENT(OUT) :: nxout
  INTEGER, INTENT(OUT) :: nyout
  INTEGER, INTENT(OUT) :: nzout
  INTEGER, INTENT(OUT) :: nzsoilout
  INTEGER, INTENT(OUT) :: nstypsout
  INTEGER, INTENT(OUT) :: istatus

!------------------------------------------------------------------------
!
!  Misc. Local variables
!
!------------------------------------------------------------------------

  INTEGER :: dimid

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begining of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Get ARPS dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_INQ_DIMID(ncid,'x_stag',dimid)
  CALL net_check_error(istatus,'NF_INQ_DIMID in netreaddims')
  istatus = NF_INQ_DIMLEN(ncid,dimid,nxout)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netreaddims')

  istatus = NF_INQ_DIMID(ncid,'y_stag',dimid)
  CALL net_check_error(istatus,'NF_INQ_DIMID in netreaddims')
  istatus = NF_INQ_DIMLEN(ncid,dimid,nyout)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netreaddims')

  istatus = NF_INQ_DIMID(ncid,'z_stag',dimid)
  CALL net_check_error(istatus,'NF_INQ_DIMID in netreaddims')
  istatus = NF_INQ_DIMLEN(ncid,dimid,nzout)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netreaddims')

  istatus = NF_INQ_DIMID(ncid,'zsoil',dimid)
  CALL net_check_error(istatus,'NF_INQ_DIMID in netreaddims')
  istatus = NF_INQ_DIMLEN(ncid,dimid,nzsoilout)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netreaddims')

  istatus = NF_INQ_DIMID(ncid,'nstyp',dimid)
  CALL net_check_error(istatus,'NF_INQ_DIMID in netreaddims')
  istatus = NF_INQ_DIMLEN(ncid,dimid,nstypsout)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netreaddims')

  RETURN
END SUBROUTINE netreaddims
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE netreadatts                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netreadatts(ncid,grdbas,runnamein,nocmntin,cmntin,dx,dy,     &
                       yearin,monthin,dayin,hourin,minutein,secondin,   &
                       thisdmpin,tstopin, &
                       mapproj,sclfct,trulat1,trulat2,trulon,latitudin, &
                       ctrlat,ctrlon,xgrdorg,ygrdorg,umovein,vmovein,   &
                       ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice,      &
                       rhosnow,rhogrpl,rhohail,alpharain,alphaice,      &
                       alphasnow,alphagrpl,alphahail,                   &
                       grdflg,basflg,varflg,mstflg,iceflg,trbflg,       &
                       sfcflg,rainflg,landflg,totflg,tkeflg,            &
                       prcflg,radflg,flxflg,snowflg,                    &
                       nscalarin,nqscalarin,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE
!
!   Retieve ARPS grib information from the NetCDF file which are stored
!   as Global attributes.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'
  INCLUDE 'globcst.inc'

  INTEGER,           INTENT(IN)  :: ncid
  INTEGER,           INTENT(IN)  :: grdbas

  CHARACTER(LEN=80), INTENT(OUT) :: runnamein
  INTEGER,           INTENT(OUT) :: nocmntin
  CHARACTER(LEN=80), INTENT(OUT) :: cmntin(50)
  REAL,              INTENT(OUT) :: dx, dy
  INTEGER,           INTENT(OUT) :: yearin, monthin, dayin, hourin,minutein,secondin
  REAL,              INTENT(OUT) :: thisdmpin, tstopin
  INTEGER,           INTENT(OUT) :: mapproj
  REAL,              INTENT(OUT) :: sclfct, trulat1, trulat2, trulon, latitudin
  REAL,              INTENT(OUT) :: ctrlat, ctrlon,  xgrdorg, ygrdorg
  REAL,              INTENT(OUT) :: umovein,  vmovein
  REAL,              INTENT(OUT) :: ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice
  REAL,              INTENT(OUT) :: rhosnow,rhogrpl,rhohail,alpharain,alphaice
  REAL,              INTENT(OUT) :: alphasnow,alphagrpl,alphahail
  INTEGER,           INTENT(OUT) :: grdflg, basflg,  varflg, mstflg,  iceflg
  INTEGER,           INTENT(OUT) :: trbflg, sfcflg, rainflg, landflg, totflg
  INTEGER,           INTENT(OUT) :: tkeflg, prcflg,  radflg, flxflg,  snowflg
  INTEGER,           INTENT(OUT) :: nscalarin
  INTEGER,           INTENT(OUT) :: nqscalarin(nscalar)
  INTEGER,           INTENT(OUT) :: istatus

  INTEGER           :: n, nq
  CHARACTER(LEN=80) :: tmpstr
  CHARACTER(LEN=1)  :: ach

  CHARACTER(LEN=17), PARAMETER :: fmtver500 = '005.10 NetCDF 3.0'
  CHARACTER(LEN=17), PARAMETER :: fmtver530 = '005.30 NetCDF 3.0'

  INTEGER,           PARAMETER :: intver500 = 500, intver530 = 530
  INTEGER :: intverin

  CHARACTER(LEN=40)  :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = NF_GET_ATT_TEXT(ncid,NF_GLOBAL,'FMTVER',tmpstr)
  CALL net_check_error(istatus,'NF_GET_ATT_TEXT (FMTVER) in netreadatts')

  IF (tmpstr(1:17) == fmtver500) THEN
    intverin = intver500
  ELSE IF (tmpstr(1:17) == fmtver530) THEN
    intverin = intver530
  ELSE
    WRITE(6,'(/1x,a,a,a/)')                                             &
              'Incoming data format, fmtverin = ',TRIM(tmpstr),         &
              ', not found. The Job stopped.'
    CALL arpsstop('arpstop called from netreadatts. ',1)
  END IF

  !
  ! Annotation
  !
  runnamein(:) = ' '
  istatus = NF_GET_ATT_TEXT(ncid,NF_GLOBAL,'RUNNAME',runnamein)
  CALL net_check_error(istatus,'NF_GET_ATT_TEXT (RUNNAME) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'nocmnt',nocmntin)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (nocmnt) in netreadatts')

  tmpstr = ' '
  DO n = 1, nocmntin
    WRITE(tmpstr,'(a,I2.2)') 'cmnt',n
    istatus = NF_GET_ATT_TEXT(ncid,NF_GLOBAL,TRIM(tmpstr),cmntin(n))
    CALL net_check_error(istatus,'NF_GET_ATT_TEXT ('//TRIM(tmpstr)//') in netreadatts')
  END DO

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DX',dx)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (DX) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DY',dy)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (DY) in netreadatts')

  !
  ! Date & time
  !
  istatus = NF_GET_ATT_TEXT(ncid,NF_GLOBAL,'INITIAL_TIME',tmpstr)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (INITIAL_TIME) in netreadatts')
  READ(tmpstr,'(I4.4,5(a,I2.2))') yearin,ach,monthin,ach,dayin,ach,     &
                                  hourin,ach,minutein,ach,secondin

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TSTOP',tstopin)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (TSTOP) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'THISDMP',thisdmpin)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (THISDMP) in netreadatts')

  !
  ! Map projection
  !
  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'MAPPROJ',mapproj)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (MAPPROJ) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'SCLFCT',sclfct)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (SCLFCT) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT1',trulat1)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (TRUELAT1) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT2',trulat2)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (TRULAT2) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELON',trulon)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (TRUELON) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'LATITUD',latitudin)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (LATITUD) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLAT',ctrlat)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (CTRLAT) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLON',ctrlon)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (CTRLON) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'XGRDORG',xgrdorg)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (XGRDORG) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'YGRDORG',ygrdorg)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (YGRDORG) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'UMOVE',umovein)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (UMOVE) in netreadatts')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'VMOVE',vmovein)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL (VMOVE) in netreadatts')

  IF (intverin < intver530) THEN
    ntcloud = 0.
    n0rain = 0.
    n0snow = 0.
    n0grpl = 0.
    n0hail = 0.
    rhoice = 0.
    rhosnow = 0.
    rhogrpl = 0.
    rhohail = 0.0
    alpharain = 0.
    alphaice = 0.
    alphasnow = 0.
    alphagrpl = 0.
    alphahail = 0.
  ELSE
    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'NTCLOUD',ntcloud)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (NTCLOUD) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'N0RAIN',n0rain)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (N0RAIN) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'N0SNOW',n0snow)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (N0SNOW) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'N0GRPL',n0grpl)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (N0GRPL) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'N0HAIL',n0hail)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (N0HAIL) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'RHOICE',rhoice)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (RHOICE) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'RHOSNOW',rhosnow)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (RHOSNOW) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'RHOGRPL',rhogrpl)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (RHOGRPL) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'RHOHAIL',rhohail)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (RHOHAIL) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'ALPHARAIN',alpharain)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (ALPHARAIN) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'ALPHAICE',alphaice)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (ALPHAICE) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'ALPHASNOW',alphasnow)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (ALPHASNOW) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'ALPHAGRPL',alphagrpl)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (ALPHAGRPL) in netreadatts')

    istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'ALPHAHAIL',alphahail)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (ALPHAHAIL) in netreadatts')

  END IF

  !
  ! Flags
  !
  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'GRDFLG',grdflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (GRDFLG) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'BASFLG',basflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (BASFLG) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'VARFLG',varflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (VARGLG) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'MSTFLG',mstflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (MSTFLG) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'ICEFLG',iceflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (ICEFLG) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'TRBFLG',trbflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (TRBFLG) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'SFCFLG',sfcflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (SFCFLG) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'RAINFLG',rainflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (RAINFLG) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'LANDFLG',landflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (LANDFLG) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'TOTFLG',totflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (TOTFLG) in netreadatts')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'TKEFLG',tkeflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT (TKEFLG) in netreadatts')

  IF (totflg == 1) THEN

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'PRCFLG',prcflg)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (PRCFLG) in netreadatts')

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'RADFLG',radflg)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (RADFLG) in netreadatts')

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'FLXFLG',flxflg)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (FLXFLG) in netreadatts')

    istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'SNOWFLG',snowflg)
    CALL net_check_error(istatus,'NF_GET_ATT_INT (SNOWFLG) in netreadatts')

  END IF

  IF (grdbas /= 1) THEN

    IF (intverin < intver530) THEN

      nqscalarin(:) = 0
      nscalarin = 0

      IF (mstflg == 1) THEN
        nscalarin = 2
        IF (P_QC > 0) nqscalarin(P_QC) = 1
        IF (P_QR > 0) nqscalarin(P_QR) = 2

        IF (iceflg == 1) THEN
          nscalarin = 5
          IF (P_QI > 0) nqscalarin(P_QI) = 3
          IF (P_QS > 0) nqscalarin(P_QS) = 4
          IF (P_QH > 0) nqscalarin(P_QH) = 5
        END IF

      END IF

    ELSE

      istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'nscalar',nscalarin)
      CALL net_check_error(istatus,'NF_GET_ATT_INT (nscalar) in netreadatts')

      tmpstr = ' '
      DO nq = 1,nscalar
        tmpstr = upcase(qnames(nq))
        istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'P_'//TRIM(tmpstr),&
                                 nqscalarin(nq))
        CALL net_check_error(istatus,'NF_GET_ATT_INT (P_'//TRIM(tmpstr)//') in netreadatts')
      END DO

    END IF

  END IF

  RETURN
END SUBROUTINE netreadatts
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE netreadattstr              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netreadattstr(ncid,attname,attstr,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE
!
!   Retieve ARPS global string from the NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER,           INTENT(IN)  :: ncid
  CHARACTER(*),      INTENT(OUT) :: attname
  CHARACTER(*),      INTENT(OUT) :: attstr
  INTEGER,           INTENT(OUT) :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = NF_GET_ATT_TEXT(ncid,NF_GLOBAL,attname,attstr)
  CALL net_check_error(istatus,'NF_GET_ATT_TEXT ('//attname//') in netreadattstr')

  RETURN
END SUBROUTINE netreadattstr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE netreadatti                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netreadatti(ncid,attname,attval,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE
!
!   Retieve ARPS integer global attributine from the NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER,           INTENT(IN)  :: ncid
  CHARACTER(*),      INTENT(OUT) :: attname
  INTEGER,           INTENT(OUT) :: attval
  INTEGER,           INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,attname,attval)
  CALL net_check_error(istatus,'NF_GET_ATT_INT ('//attname//') in netreadatti')

  RETURN
END SUBROUTINE netreadatti

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE netreadattr                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netreadattr(ncid,attname,attval,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE
!
!   Retieve ARPS integer global attributine from the NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER,           INTENT(IN)  :: ncid
  CHARACTER(*),      INTENT(OUT) :: attname
  REAL,              INTENT(OUT) :: attval
  INTEGER,           INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,attname,attval)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL ('//attname//') in netreadatti')

  RETURN
END SUBROUTINE netreadattr
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE netread1d                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netread1d(ncid,packed,itime,varname,nx,var1d)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 1D array from the ARPS NetCDF file.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: packed
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: nx
  REAL,             INTENT(OUT) :: var1d(nx)

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER           :: istatus
  INTEGER           :: varid
  CHARACTER(LEN=20) :: namein
  INTEGER           :: vartype, ndims,natts,dimlen
  INTEGER           :: dimids(5)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(6,FMT='(1x,2a)',ADVANCE='NO')     &
                        '  Reading 1D NetCDF variable ', varname
  !
  ! Get variable ID
  !
  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in netread1d')

  !
  ! Do some checks
  !
  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL net_check_error(istatus,'NF_INQ_VAR in netread1d')

  IF(vartype /= NF_FLOAT) THEN      ! Data type
    WRITE(6,'(1x,3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF
                                    ! Data rank
  IF((ndims /= 2 .AND. itime > 0) .OR. (ndims /= 1 .AND. itime == 0) ) THEN
    WRITE(6,'(1x,3a)') 'Variable ', varname, ' is not a 1D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF
                                    ! X dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread1d')
  IF(dimlen /= nx) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'First dimension of variable ', varname, &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF

  IF (itime > 0) THEN               ! Record No. if applied
    istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlen)
    CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread1d')
    IF(dimlen < itime) THEN
      WRITE(6,'(1x,a,I3,a,I3)') 'The total records number is ', dimlen, &
                      ' however, the required time level is ',itime
      STOP 'itime_tool_large'
    END IF
  END IF

  IF (itime > 0) THEN
    istatus = NF_GET_VARA_REAL(ncid,varid,(/1,itime/),(/nx,1/),var1d)
    CALL net_check_error(istatus,'NF_GET_VARA_REAL in netread1d')
  ELSE
    istatus = NF_GET_VARA_REAL(ncid,varid,(/1/),(/nx/),var1d)
    CALL net_check_error(istatus,'NF_GET_VARA_REAL in netread1d')
  END IF

  WRITE(6,'(a)') '         === DONE ==='

  RETURN
END SUBROUTINE netread1d
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE netread2d                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netread2d(ncid,packed,itime,varname,nx,ny,var2d)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 2D array from the ARPS NetCDF file.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: packed
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: nx
  INTEGER,          INTENT(IN)  :: ny
  REAL,             INTENT(OUT) :: var2d(nx,ny)

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER           :: istatus
  INTEGER           :: varid
  CHARACTER(LEN=20) :: namein
  INTEGER           :: vartype, ndims,natts,dimlen
  INTEGER           :: dimids(5)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(6,FMT='(1x,2a)',ADVANCE='NO')     &
                        '  Reading 2D NetCDF variable ', varname
  !
  ! Get variable ID
  !
  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in netread2d')

  !
  ! Do some checks
  !
  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL net_check_error(istatus,'NF_INQ_VAR in netread2d')

  IF(vartype /= NF_FLOAT) THEN      ! Data type
    WRITE(6,'(1x,3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF
                                    ! Data rank
  IF((ndims /= 3 .AND. itime > 0) .OR. (ndims /= 2 .AND. itime == 0) ) THEN
    WRITE(6,'(1x,3a)') 'Variable ', varname, ' is not a 2D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF
                                    ! X dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread2d')
  IF(dimlen /= nx) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'First dimension of variable ', varname, &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF
                                    ! Y dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread2d')
  IF(dimlen /= ny) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'Second dimension of variable ',varname, &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF

  IF (itime > 0) THEN               ! Record No. if applied
    istatus = NF_INQ_DIMLEN(ncid,dimids(3),dimlen)
    CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread2d')
    IF(dimlen < itime) THEN
      WRITE(6,'(1x,a,I3,a,I3)') 'The total records number is ', dimlen, &
                      ' however, the required time level is ',itime
      STOP 'itime_tool_large'
    END IF
  END IF

  IF (itime > 0) THEN
    istatus = NF_GET_VARA_REAL(ncid,varid,(/1,1,itime/),(/nx,ny,1/),var2d)
    CALL net_check_error(istatus,'NF_GET_VARA_REAL in netread2d')
  ELSE
    istatus = NF_GET_VARA_REAL(ncid,varid,(/1,1/),(/nx,ny/),var2d)
    CALL net_check_error(istatus,'NF_GET_VARA_REAL in netread2d')
  END IF

  WRITE(6,'(a)') '         === DONE ==='

  RETURN
END SUBROUTINE netread2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE netread2di                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netread2di(ncid,packed,itime,varname,nx,ny,var2d)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 2D INTEGER array from the ARPS NetCDF file.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: packed
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: nx
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(OUT) :: var2d(nx,ny)

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER           :: istatus
  INTEGER           :: varid
  CHARACTER(LEN=20) :: namein
  INTEGER           :: vartype, ndims,natts,dimlen
  INTEGER           :: dimids(5)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(6,FMT='(1x,2a)',ADVANCE='NO')     &
                      '  Reading 2D integer NetCDF variable ', varname
  !
  ! Get variable ID
  !
  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in netread2di')

  !
  ! Do some checks
  !
  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL net_check_error(istatus,'NF_INQ_VAR in netread2di')

  IF(vartype /= NF_INT) THEN        ! Data type
    WRITE(6,'(1x,3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF
                                    ! Data rank
  IF((ndims /= 3 .AND. itime > 0) .OR. (ndims /= 2 .AND. itime == 0) ) THEN
    WRITE(6,'(1x,3a)') 'Variable ', varname, ' is not a 2D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF
                                    ! X dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread2di')
  IF(dimlen /= nx) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'First dimension of variable ', varname, &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF
                                    ! Y dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread2di')
  IF(dimlen /= ny) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'Second dimension of variable ',varname, &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF

  IF (itime > 0) THEN               ! Record No. if applied
    istatus = NF_INQ_DIMLEN(ncid,dimids(3),dimlen)
    CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread2di')
    IF(dimlen < itime) THEN
      WRITE(6,'(1x,a,I3,a,I3)') 'The total records number is ', dimlen, &
                      ' however, the required time level is ',itime
      STOP 'itime_tool_large'
    END IF
  END IF

  IF (itime > 0) THEN
    istatus = NF_GET_VARA_INT(ncid,varid,(/1,1,itime/),(/nx,ny,1/),var2d)
    CALL net_check_error(istatus,'NF_GET_VARA_REAL in netread2di')
  ELSE
    istatus = NF_GET_VARA_INT(ncid,varid,(/1,1/),(/nx,ny/),var2d)
    CALL net_check_error(istatus,'NF_GET_VARA_REAL in netread2di')
  END IF

  WRITE(6,'(a)') '         === DONE ==='

  RETURN
END SUBROUTINE netread2di
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE netread3d                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netread3d(ncid,packed,itime,varname,nx,ny,nz,var3d)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 3D array from the ARPS NetCDF file.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: packed
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: nx
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(IN)  :: nz
  REAL,             INTENT(OUT) :: var3d(nx,ny,nz)

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER           :: istatus
  INTEGER           :: varid
  CHARACTER(LEN=20) :: namein
  INTEGER           :: vartype, ndims,natts,dimlen
  INTEGER           :: dimids(5)

  CHARACTER(LEN=20) :: fmtstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(6,FMT='(1x,2a)',ADVANCE='NO')     &
                        '  Reading 3D NetCDF variable ', TRIM(varname)
  !
  ! Get variable ID
  !
  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in netread3d')

  !
  ! Do some checks
  !
  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL net_check_error(istatus,'NF_INQ_VAR in netread3d')

  IF(vartype /= NF_FLOAT) THEN      ! Data type
    WRITE(6,'(1x,3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF
                                    ! Data rank
  IF((ndims /= 4 .AND. itime > 0) .OR. (ndims /= 3 .AND. itime == 0) ) THEN
    WRITE(6,'(1x,3a)') 'Variable ', varname, ' is not a 2D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF
                                    ! X dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread3d')
  IF(dimlen /= nx) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'First dimension of variable ', varname, &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF
                                    ! Y dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread3d')
  IF(dimlen /= ny) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'Second dimension of variable ',varname, &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF
                                    ! Z dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(3),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread3d')
  IF(dimlen /= nz) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'Third dimension of variable ',varname, &
                    ' is ',dimlen, ' and it should be ',nz
    STOP 'WRONG_DIM_length'
  END IF

  IF (itime > 0) THEN               ! Record No. if applied
    istatus = NF_INQ_DIMLEN(ncid,dimids(4),dimlen)
    CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread3d')
    IF(dimlen < itime) THEN
      WRITE(6,'(1x,a,I3,a,I3)') 'The total records number is ', dimlen, &
                      ' however, the required time level is ',itime
      STOP 'itime_tool_large'
    END IF
  END IF

  IF (itime > 0) THEN
    istatus = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,itime/),              &
                               (/nx,ny,nz,1/),var3d)
  ELSE
    istatus = NF_GET_VARA_REAL(ncid,varid,(/1,1,1/),                    &
                               (/nx,ny,nz/),var3d)
  END IF
  CALL net_check_error(istatus,'NF_GET_VARA_REAL in netread3d')

  WRITE(fmtstr,'(a,I2.2,a)') '(',15-LEN_TRIM(varname),'x,a)'
  WRITE(6,FMT=fmtstr) ' === DONE ==='

  RETURN
END SUBROUTINE netread3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE netread3di                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netread3di(ncid,packed,itime,varname,nx,ny,nz,var3d)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 3D integer array from the ARPS NetCDF file.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: packed
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: nx
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(IN)  :: nz
  INTEGER,          INTENT(OUT) :: var3d(nx,ny,nz)

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER           :: istatus
  INTEGER           :: varid
  CHARACTER(LEN=20) :: namein
  INTEGER           :: vartype, ndims,natts,dimlen
  INTEGER           :: dimids(5)

  CHARACTER(LEN=20) :: fmtstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(6,FMT='(1x,2a)',ADVANCE='NO')     &
                        '  Reading 3D integer NetCDF variable ', varname
  !
  ! Get variable ID
  !
  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in netread3di')

  !
  ! Do some checks
  !
  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL net_check_error(istatus,'NF_INQ_VAR in netread3di')

  IF(vartype /= NF_INT) THEN        ! Data type
    WRITE(6,'(1x,3a)') 'Variable ',varname, ' is not INTEGER.'
    STOP 'WRONG_VAR_TYPE'
  END IF
                                    ! Data rank
  IF((ndims /= 4 .AND. itime > 0) .OR. (ndims /= 3 .AND. itime == 0) ) THEN
    WRITE(6,'(1x,3a)') 'Variable ', varname, ' is not a 2D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF
                                    ! X dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread3di')
  IF(dimlen /= nx) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'First dimension of variable ', varname, &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF
                                    ! Y dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread3di')
  IF(dimlen /= ny) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'Second dimension of variable ',varname, &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF
                                    ! Z dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(3),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread3di')
  IF(dimlen /= nz) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'Third dimension of variable ',varname, &
                    ' is ',dimlen, ' and it should be ',nz
    STOP 'WRONG_DIM_length'
  END IF

  IF (itime > 0) THEN               ! Record No. if applied
    istatus = NF_INQ_DIMLEN(ncid,dimids(4),dimlen)
    CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread3di')
    IF(dimlen < itime) THEN
      WRITE(6,'(1x,a,I3,a,I3)') 'The total records number is ', dimlen, &
                      ' however, the required time level is ',itime
      STOP 'itime_tool_large'
    END IF
  END IF

  IF (itime > 0) THEN
    istatus = NF_GET_VARA_INT(ncid,varid,(/1,1,1,itime/),              &
                               (/nx,ny,nz,1/),var3d)
  ELSE
    istatus = NF_GET_VARA_INT(ncid,varid,(/1,1,1/),                    &
                               (/nx,ny,nz/),var3d)
  END IF
  CALL net_check_error(istatus,'NF_GET_VARA_REAL in netread3di')

  WRITE(fmtstr,'(a,I2.2,a)') '(',14-LEN_TRIM(varname),'x,a)'
  WRITE(6,FMT=fmtstr) '         === DONE ==='

  RETURN
END SUBROUTINE netread3di
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE netread4d                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netread4d(ncid,packed,itime,varname,nx,ny,nz,nn,var4d)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 3D array from the ARPS NetCDF file.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: packed
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: nx
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(IN)  :: nz
  INTEGER,          INTENT(IN)  :: nn
  REAL,             INTENT(OUT) :: var4d(nx,ny,nz,nn)

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER           :: istatus
  INTEGER           :: varid
  CHARACTER(LEN=20) :: namein
  INTEGER           :: vartype, ndims,natts,dimlen
  INTEGER           :: dimids(5)

  CHARACTER(LEN=20) :: fmtstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(6,FMT='(1x,2a)',ADVANCE='NO')     &
                        '  Reading 4D NetCDF variable ', varname
  !
  ! Get variable ID
  !
  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in netread4d')

  !
  ! Do some checks
  !
  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL net_check_error(istatus,'NF_INQ_VAR in netread4d')

  IF(vartype /= NF_FLOAT) THEN      ! Data type
    WRITE(6,'(1x,3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF
                                    ! Data rank
  IF((ndims /= 5 .AND. itime > 0) .OR. (ndims /= 4 .AND. itime == 0) ) THEN
    WRITE(6,'(1x,3a)') 'Variable ', varname, ' is not a 2D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF
                                    ! X dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread4d')
  IF(dimlen /= nx) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'First dimension of variable ', varname, &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF
                                    ! Y dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread4d')
  IF(dimlen /= ny) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'Second dimension of variable ',varname, &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF
                                    ! Z dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(3),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread4d')
  IF(dimlen /= nz) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'Third dimension of variable ',varname, &
                    ' is ',dimlen, ' and it should be ',nz
    STOP 'WRONG_DIM_length'
  END IF
                                    ! nstyps dimension length
  istatus = NF_INQ_DIMLEN(ncid,dimids(4),dimlen)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread4d')
  IF(dimlen /= nn) THEN
    WRITE(6,'(1x,3a,I3,a,I3)') 'Fourth dimension of variable ',varname, &
                    ' is ',dimlen, ' and it should be ',nn
    STOP 'WRONG_DIM_length'
  END IF

  IF (itime > 0) THEN               ! Record No. if applied
    istatus = NF_INQ_DIMLEN(ncid,dimids(5),dimlen)
    CALL net_check_error(istatus,'NF_INQ_DIMLEN in netread4d')
    IF(dimlen < itime) THEN
      WRITE(6,'(1x,a,I3,a,I3)') 'The total records number is ', dimlen, &
                      ' however, the required time level is ',itime
      STOP 'itime_tool_large'
    END IF
  END IF

  IF (itime > 0) THEN
    istatus = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,1,itime/),              &
                               (/nx,ny,nz,nn,1/),var4d)
  ELSE
    istatus = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,1/),                    &
                               (/nx,ny,nz,nn/),var4d)
  END IF
  CALL net_check_error(istatus,'NF_GET_VARA_REAL in netread4d')

  WRITE(fmtstr,'(a,I2.2,a)') '(',15-LEN_TRIM(varname),'x,a)'
  WRITE(6,FMT=fmtstr) '         === DONE ==='

  RETURN
END SUBROUTINE netread4d
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE netreadTime                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE netreadTime(ncid,itime,varname,time)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in current Valid time from the ARPS NetCDF file.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  REAL,             INTENT(OUT) :: time

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER           :: istatus
  INTEGER           :: varid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(6,FMT='(1x,2a)',ADVANCE='NO')     &
                        '  Reading data valid time ', varname
  !
  ! Get variable ID
  !
  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL net_check_error(istatus,'NF_INQ_VARID in netreadTime')

  IF (itime > 0) THEN
    istatus = NF_GET_VARA_REAL(ncid,varid,(/itime/),(/1/),time)
  ELSE
    WRITE(6,'(1x,a)') 'There is not CURTIME variable in the data file.'
    STOP
  END IF
  CALL net_check_error(istatus,'NF_GET_VARA_REAL in netreadTime')

  WRITE(6,'(a)') '         === DONE ==='

  RETURN
END SUBROUTINE netreadTime
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_define_trn                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_define_trn(ncid,nx,ny,dx,dy,mapproj,sclfct,              &
                          trulat1,trulat2,trulon,ctrlat,ctrlon,         &
                          istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Define ARPS terrain data file attributes and variables. After this call
!     The netCDF file should be in DATA mode.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (08/10/2004)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  USE arps_netio_metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(IN)  :: nx, ny
  REAL,    INTENT(IN)  :: dx, dy
  INTEGER, INTENT(IN)  :: mapproj
  REAL,    INTENT(IN)  :: sclfct
  REAL,    INTENT(IN)  :: trulat1, trulat2, trulon
  REAL,    INTENT(IN)  :: ctrlat,  ctrlon

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER           :: varid
  INTEGER           :: dimx_id,dimy_id

  INTEGER           :: oldfillmode

  CHARACTER(LEN=80) :: tmpstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Define dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_DEF_DIM(ncid,'x',nx-1,dimx_id)
  CALL net_check_error(istatus,'net_define_trn')

  istatus = NF_DEF_DIM(ncid,'y',ny-1,dimy_id)
  CALL net_check_error(istatus,'net_define_trn')

!-----------------------------------------------------------------------
!
! Define global attributes
!
!-----------------------------------------------------------------------

  tmpstr   = 'ARPS 5.2 terrain data (ARPSTRN)'

  CALL netwrt_general_att(ncid,0,tmpstr,nx,ny,dx,dy,mapproj,sclfct,     &
                          trulat1,trulat2,trulon,ctrlat,ctrlon,istatus)

  ! do not fill, will set values explicitly later. Improve performance

  istatus = NF_SET_FILL(ncid,NF_NOFILL,oldfillmode)

!-----------------------------------------------------------------------
!
! Define variable arrays
!
!-----------------------------------------------------------------------

  istatus = NF_DEF_VAR(ncid,'HTERAIN',NF_FLOAT,2,(/dimx_id,dimy_id/),varid)
  CALL net_define_var_meta(ncid,varid,'REAL',trnmeta%hterrain)

!-----------------------------------------------------------------------
!
! End NetCDF file DEFINE mode
!
!-----------------------------------------------------------------------

  istatus = NF_ENDDEF(ncid)
  CALL net_check_error(istatus,'net_define_trn')

  RETURN
END SUBROUTINE net_define_trn
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_get_trn                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_get_trn(ncid,nx,ny,dx,dy,mapproj,sclfct,              &
                        trulat1,trulat2,trulon,ctrlat,ctrlon,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Extract NetCDF file dimensions and attributes from ARPS terrain
!     data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (08/11/2004)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(OUT) :: nx, ny
  REAL,    INTENT(OUT) :: dx, dy
  INTEGER, INTENT(OUT) :: mapproj
  REAL,    INTENT(OUT) :: sclfct
  REAL,    INTENT(OUT) :: trulat1, trulat2, trulon
  REAL,    INTENT(OUT) :: ctrlat,  ctrlon

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER           :: varid
  INTEGER           :: dimx_id,dimy_id

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Get dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_INQ_DIMID(ncid,'x',dimx_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_trn')
  istatus = NF_INQ_DIMLEN(ncid,dimx_id,nx)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_trn')

  istatus = NF_INQ_DIMID(ncid,'y',dimy_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_trn')
  istatus = NF_INQ_DIMLEN(ncid,dimy_id,ny)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_trn')

  nx = nx + 1
  ny = ny + 1
!-----------------------------------------------------------------------
!
! Get global attributes
!
!-----------------------------------------------------------------------

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DX',dx)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_trn')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DY',dy)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_trn')

  !
  ! Map projection
  !
  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'MAPPROJ',mapproj)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_trn')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'SCLFCT',sclfct)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_trn')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT1',trulat1)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_trn')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT2',trulat2)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_trn')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELON',trulon)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_trn')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLAT',ctrlat)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_trn')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLON',ctrlon)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_trn')

  RETURN
END SUBROUTINE net_get_trn
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_define_sfc                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_define_sfc(ncid,nx,ny,nstyps,dx,dy,mapproj,sclfct,       &
               trulat1,trulat2,trulon,ctrlat,ctrlon,                    &
               stypflg,vtypflg,laiflg,rfnsflg,vegflg,ndviflg,           &
               istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Define ARPS surface data file attributes and variables. After this call
!     The netCDF file should be in DATA mode.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (08/12/2004)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  use arps_netio_metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(IN)  :: nx, ny, nstyps
  REAL,    INTENT(IN)  :: dx, dy
  INTEGER, INTENT(IN)  :: mapproj
  REAL,    INTENT(IN)  :: sclfct
  REAL,    INTENT(IN)  :: trulat1, trulat2, trulon
  REAL,    INTENT(IN)  :: ctrlat,  ctrlon
  INTEGER, INTENT(IN)  :: stypflg, vtypflg, laiflg
  INTEGER, INTENT(IN)  :: rfnsflg, vegflg,  ndviflg

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER           :: varid
  INTEGER           :: dimx_id, dimy_id, dimn_id

  INTEGER           :: oldfillmode

  CHARACTER(LEN=80) :: tmpstr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
! Define dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_DEF_DIM(ncid,'x',nx-1,dimx_id)
  CALL net_check_error(istatus,'net_define_sfc')

  istatus = NF_DEF_DIM(ncid,'y',ny-1,dimy_id)
  CALL net_check_error(istatus,'net_define_sfc')

  istatus = NF_DEF_DIM(ncid,'nstyp',nstyps,dimn_id)
  CALL net_check_error(istatus,'net_define_sfc')

!-----------------------------------------------------------------------
!
! Define global attributes
!
!-----------------------------------------------------------------------

  tmpstr   = 'ARPS 5.2 surface characteristics data (ARPSSFC)'

  CALL netwrt_general_att(ncid,0,tmpstr,nx,ny,dx,dy,mapproj,sclfct,     &
                          trulat1,trulat2,trulon,ctrlat,ctrlon,istatus)
  !
  ! Flags
  !
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'STYPFLG',NF_INT,1,stypflg)
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'VTYPFLG',NF_INT,1,vtypflg)
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'LAIFLG', NF_INT,1,laiflg)
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'RFNSFLG',NF_INT,1,rfnsflg)
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'VEGFLG', NF_INT,1,vegflg)
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NDVIFLG',NF_INT,1,ndviflg)

  ! do not fill, will set values explicitly later. Improve performance

  istatus = NF_SET_FILL(ncid,NF_NOFILL,oldfillmode)

!-----------------------------------------------------------------------
!
! Define variable arrays
!
!-----------------------------------------------------------------------

  IF (stypflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'SOILTYP',NF_INT,3,                       &
                                     (/dimx_id,dimy_id,dimn_id/),varid)
    CALL net_define_var_meta(ncid,varid,'INT',sfcmeta%soiltyp)

    istatus = NF_DEF_VAR(ncid,'STYPFRCT',NF_FLOAT,3,                    &
                                     (/dimx_id,dimy_id,dimn_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',sfcmeta%stypfrct)

  END IF

  IF (vtypflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'VEGTYP',NF_INT,2,(/dimx_id,dimy_id/),varid)
    CALL net_define_var_meta(ncid,varid,'INT',sfcmeta%vegtyp)
  END IF

  IF (laiflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'LAI',NF_FLOAT,2,(/dimx_id,dimy_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',sfcmeta%lai)
  END IF

  IF (rfnsflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'ROUFNS',NF_FLOAT,2,(/dimx_id,dimy_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',sfcmeta%roufns)
  END IF

  IF (vegflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'VEG',NF_FLOAT,2,(/dimx_id,dimy_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',sfcmeta%veg)
  END IF

  IF (ndviflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'NDVI',NF_FLOAT,2,(/dimx_id,dimy_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',sfcmeta%ndvi)
  END IF

!-----------------------------------------------------------------------
!
! End NetCDF file DEFINE mode
!
!-----------------------------------------------------------------------

  istatus = NF_ENDDEF(ncid)
  CALL net_check_error(istatus,'net_define_sfc')

  RETURN
END SUBROUTINE net_define_sfc
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_get_sfc                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_get_sfc(ncid,nx,ny,nstyps,dx,dy,                         &
               mapproj,sclfct,trulat1,trulat2,trulon,ctrlat,ctrlon,     &
               stypflg,vtypflg,laiflg,rfnsflg,vegflg,ndviflg,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Extract NetCDF file dimensions and attributes from ARPS surface
!     data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (08/12/2004)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(OUT) :: nx, ny, nstyps
  REAL,    INTENT(OUT) :: dx, dy
  INTEGER, INTENT(OUT) :: mapproj
  REAL,    INTENT(OUT) :: sclfct
  REAL,    INTENT(OUT) :: trulat1, trulat2, trulon
  REAL,    INTENT(OUT) :: ctrlat,  ctrlon
  INTEGER, INTENT(OUT) :: stypflg, vtypflg,  laiflg
  INTEGER, INTENT(OUT) :: rfnsflg, vegflg,   ndviflg

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dimx_id,dimy_id,dimn_id

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Get dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_INQ_DIMID(ncid,'x',dimx_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_sfc')
  istatus = NF_INQ_DIMLEN(ncid,dimx_id,nx)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_sfc')

  istatus = NF_INQ_DIMID(ncid,'y',dimy_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_sfc')
  istatus = NF_INQ_DIMLEN(ncid,dimy_id,ny)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_sfc')

  istatus = NF_INQ_DIMID(ncid,'nstyp',dimn_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_sfc')
  istatus = NF_INQ_DIMLEN(ncid,dimn_id,nstyps)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_sfc')

  nx = nx + 1
  ny = ny + 1
!-----------------------------------------------------------------------
!
! Get global attributes
!
!-----------------------------------------------------------------------

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DX',dx)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_sfc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DY',dy)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_sfc')

  !
  ! Map projection
  !
  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'MAPPROJ',mapproj)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_sfc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'SCLFCT',sclfct)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_sfc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT1',trulat1)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_sfc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT2',trulat2)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_sfc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELON',trulon)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_sfc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLAT',ctrlat)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_sfc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLON',ctrlon)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_sfc')

  !
  ! Flags
  !
  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'STYPFLG',stypflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_sfc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'VTYPFLG',vtypflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_sfc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'LAIFLG',laiflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_sfc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'RFNSFLG',rfnsflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_sfc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'VEGFLG',vegflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_sfc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'NDVIFLG',ndviflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_sfc')

  RETURN
END SUBROUTINE net_get_sfc
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_define_soil                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_define_soil(ncid,nx,ny,nzsoil,nstyps,dx,dy,mapproj,      &
                  sclfct,trulat1,trulat2,trulon,ctrlat,ctrlon,          &
                  zpsoilflg,tsoilflg,qsoilflg,wcanpflg,                 &
                  snowdflg,stypflg,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Define ARPS soil data file attributes and variables. After this call
!     The netCDF file should be in DATA mode.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (08/13/2004)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  USE arps_netio_metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(IN)  :: nx, ny, nzsoil, nstyps
  REAL,    INTENT(IN)  :: dx, dy
  INTEGER, INTENT(IN)  :: mapproj
  REAL,    INTENT(IN)  :: sclfct
  REAL,    INTENT(IN)  :: trulat1, trulat2, trulon
  REAL,    INTENT(IN)  :: ctrlat,  ctrlon
  INTEGER, INTENT(IN)  :: zpsoilflg, tsoilflg, qsoilflg
  INTEGER, INTENT(IN)  :: wcanpflg,  snowdflg, stypflg

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER           :: varid
  INTEGER           :: dimx_id, dimy_id, dimz_id, dimn_id
  INTEGER           :: dims_id

  INTEGER           :: oldfillmode

  CHARACTER(LEN=80) :: tmpstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Define dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_DEF_DIM(ncid,'x',nx-1,dimx_id)
  CALL net_check_error(istatus,'net_define_soil')

  istatus = NF_DEF_DIM(ncid,'y',ny-1,dimy_id)
  CALL net_check_error(istatus,'net_define_soil')

  istatus = NF_DEF_DIM(ncid,'zsoil',nzsoil,dimz_id)
  CALL net_check_error(istatus,'net_define_soil')

  istatus = NF_DEF_DIM(ncid,'nstyp',nstyps,dimn_id)
  CALL net_check_error(istatus,'net_define_soil')

  istatus = NF_DEF_DIM(ncid,'nstyp_total',nstyps+1,dims_id)
  CALL net_check_error(istatus,'net_define_soil')

!-----------------------------------------------------------------------
!
! Define global attributes
!
!-----------------------------------------------------------------------

  tmpstr   = 'ARPS 5.2 Soil data'

  CALL netwrt_general_att(ncid,0,tmpstr,nx,ny,dx,dy,mapproj,sclfct,     &
                          trulat1,trulat2,trulon,ctrlat,ctrlon,istatus)
  !
  ! Flags
  !
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZPSOILFLG',NF_INT,1,zpsoilflg)
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'TSOILFLG', NF_INT,1,tsoilflg)
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QSOILFLG', NF_INT,1,qsoilflg)
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'WCANPFLG', NF_INT,1,wcanpflg)
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'SNOWDFLG', NF_INT,1,snowdflg)
  istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'STYPFLG',  NF_INT,1,stypflg)

  ! do not fill, will set values explicitly later. Improve performance

  istatus = NF_SET_FILL(ncid,NF_NOFILL,oldfillmode)

!-----------------------------------------------------------------------
!
! Define variable arrays
!
!-----------------------------------------------------------------------

  IF (zpsoilflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'ZPSOIL',NF_FLOAT,3,                      &
                         (/dimx_id,dimy_id,dimz_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',soilmeta%zpsoil)
  END IF

  IF (tsoilflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'TSOIL',NF_FLOAT,4,                       &
                         (/dimx_id,dimy_id,dimz_id,dims_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',soilmeta%tsoil)
  END IF

  IF (qsoilflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'QSOIL',NF_FLOAT,4,                       &
                         (/dimx_id,dimy_id,dimz_id,dims_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',soilmeta%qsoil)
  END IF

  IF (wcanpflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'WETCANP',NF_FLOAT,3,                     &
                         (/dimx_id,dimy_id,dims_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',soilmeta%wetcanp)
  END IF

  IF (snowdflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'SNOWDPTH',NF_FLOAT,2,(/dimx_id,dimy_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',soilmeta%snowdpth)
  END IF

  IF (stypflg /= 0) THEN
    istatus = NF_DEF_VAR(ncid,'SOILTYP',NF_INT,3,                       &
                                     (/dimx_id,dimy_id,dimn_id/),varid)
    CALL net_define_var_meta(ncid,varid,'INT',soilmeta%soiltyp)
  END IF

!-----------------------------------------------------------------------
!
! End NetCDF file DEFINE mode
!
!-----------------------------------------------------------------------

  istatus = NF_ENDDEF(ncid)
  CALL net_check_error(istatus,'net_define_soil')

  RETURN
END SUBROUTINE net_define_soil
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_get_soil                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_get_soil(ncid,nx,ny,nzsoil,nstyps,dx,dy,                 &
               mapproj,sclfct,trulat1,trulat2,trulon,ctrlat,ctrlon,     &
               zpsoilflg,tsoilflg,qsoilflg,wcanpflg,snowdflg,stypflg,   &
               istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Extract NetCDF file dimensions and attributes from ARPS soil
!     data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (08/13/2004)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(OUT) :: nx, ny, nzsoil, nstyps
  REAL,    INTENT(OUT) :: dx, dy
  INTEGER, INTENT(OUT) :: mapproj
  REAL,    INTENT(OUT) :: sclfct
  REAL,    INTENT(OUT) :: trulat1, trulat2, trulon
  REAL,    INTENT(OUT) :: ctrlat,  ctrlon
  INTEGER, INTENT(OUT) :: zpsoilflg, tsoilflg,  qsoilflg
  INTEGER, INTENT(OUT) :: wcanpflg,  snowdflg,  stypflg

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dimx_id,dimy_id,dimz_id,dimn_id,dims_id
  INTEGER :: ntotal

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Get dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_INQ_DIMID(ncid,'x',dimx_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_soil')
  istatus = NF_INQ_DIMLEN(ncid,dimx_id,nx)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_soil')

  istatus = NF_INQ_DIMID(ncid,'y',dimy_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_soil')
  istatus = NF_INQ_DIMLEN(ncid,dimy_id,ny)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_soil')

  istatus = NF_INQ_DIMID(ncid,'zsoil',dimz_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_soil')
  istatus = NF_INQ_DIMLEN(ncid,dimz_id,nzsoil)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_soil')

  istatus = NF_INQ_DIMID(ncid,'nstyp',dimn_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_soil')
  istatus = NF_INQ_DIMLEN(ncid,dimn_id,nstyps)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_soil')

  istatus = NF_INQ_DIMID(ncid,'nstyp_total',dims_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_soil')
  istatus = NF_INQ_DIMLEN(ncid,dims_id,ntotal)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_soil')

  nx = nx + 1
  ny = ny + 1
!-----------------------------------------------------------------------
!
! Get global attributes
!
!-----------------------------------------------------------------------

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DX',dx)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_soil')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DY',dy)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_soil')

  !
  ! Map projection
  !
  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'MAPPROJ',mapproj)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_soil')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'SCLFCT',sclfct)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_soil')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT1',trulat1)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_soil')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT2',trulat2)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_soil')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELON',trulon)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_soil')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLAT',ctrlat)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_soil')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLON',ctrlon)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL in net_get_soil')

  !
  ! Flags
  !
  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'ZPSOILFLG',zpsoilflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_soil')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'TSOILFLG',tsoilflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_soil')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'QSOILFLG',qsoilflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_soil')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'WCANPFLG',wcanpflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_soil')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'SNOWDFLG',snowdflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_soil')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'STYPFLG',stypflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT in net_get_soil')

  RETURN
END SUBROUTINE net_get_soil
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_define_exbc                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_define_exbc(ncid,nx,ny,nz,itime,dx,dy,dz,dzmin,strhopt,  &
             zrefsfc,dlayer1,dlayer2,zflat,strhtune,mapproj,sclfct,     &
             trulat1,trulat2,trulon,ctrlat,ctrlon,                      &
             ubcflg,vbcflg,wbcflg,ptbcflg,prbcflg,qvbcflg,              &
             nscalarout,qscalarbcflg,ctime,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Define ARPS boundary data file attributes and variables. After this call
!     The netCDF file should be in DATA mode.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (08/18/2004)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  USE arps_netio_metadata

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'
  INCLUDE 'globcst.inc'

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(IN)  :: nx, ny, nz, itime
  REAL,    INTENT(IN)  :: dx, dy, dz
  REAL,    INTENT(IN)  :: dzmin, zrefsfc, dlayer1, dlayer2, zflat
  REAL,    INTENT(IN)  :: strhtune
  INTEGER, INTENT(IN)  :: strhopt
  INTEGER, INTENT(IN)  :: mapproj
  REAL,    INTENT(IN)  :: sclfct
  REAL,    INTENT(IN)  :: trulat1, trulat2, trulon
  REAL,    INTENT(IN)  :: ctrlat,  ctrlon
  INTEGER, INTENT(IN)  :: ubcflg,  vbcflg,  wbcflg,  ptbcflg, prbcflg, qvbcflg
  INTEGER, INTENT(IN)  :: nscalarout
  INTEGER, INTENT(IN)  :: qscalarbcflg(nscalar)

  CHARACTER(LEN=*), INTENT(IN) :: ctime

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER           :: varid
  INTEGER           :: dimx_id,  dimy_id,  dimz_id, dimt_id, dims_id
  INTEGER           :: dimxs_id, dimys_id, dimzs_id

  INTEGER           :: oldfillmode

  CHARACTER(LEN=80) :: tmpstr

  INTEGER,PARAMETER :: ctime_len = 15

  INTEGER           :: nq

  CHARACTER(LEN=4) :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (itime == 1) THEN

!-----------------------------------------------------------------------
!
! Define dimensions
!
!-----------------------------------------------------------------------

    istatus = NF_DEF_DIM(ncid,'Time',NF_UNLIMITED,dimt_id)
    CALL net_check_error(istatus,'net_define_exbc')

    istatus = NF_DEF_DIM(ncid,'x_stag',nx,dimx_id)
    CALL net_check_error(istatus,'net_define_exbc')

    istatus = NF_DEF_DIM(ncid,'y_stag',ny,dimy_id)
    CALL net_check_error(istatus,'net_define_exbc')

    istatus = NF_DEF_DIM(ncid,'z_stag',nz,dimz_id)
    CALL net_check_error(istatus,'net_define_exbc')

    istatus = NF_DEF_DIM(ncid,'x',nx-1,dimxs_id)
    CALL net_check_error(istatus,'net_define_exbc')

    istatus = NF_DEF_DIM(ncid,'y',ny-1,dimys_id)
    CALL net_check_error(istatus,'net_define_exbc')

    istatus = NF_DEF_DIM(ncid,'z',nz-1,dimzs_id)
    CALL net_check_error(istatus,'net_define_exbc')

    istatus = NF_DEF_DIM(ncid,'CtimeStrLen',ctime_len,dims_id)
    CALL net_check_error(istatus,'net_define_exbc')

!-----------------------------------------------------------------------
!
! Define global attributes
!
!-----------------------------------------------------------------------

    tmpstr   = 'ARPS 5.3 EXBC data'

    CALL netwrt_general_att(ncid,0,tmpstr,nx,ny,dx,dy,mapproj,sclfct,   &
                            trulat1,trulat2,trulon,ctrlat,ctrlon,istatus)

    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'DZ',NF_FLOAT,1,dz)

    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'DZMIN',  NF_FLOAT,1,dzmin)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'ZREFSFC',NF_FLOAT,1,zrefsfc)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'DLAYER1',NF_FLOAT,1,dlayer1)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'DLAYER2',NF_FLOAT,1,dlayer2)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'ZFLAT',  NF_FLOAT,1,zflat)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'STRHTUNE',NF_FLOAT,1,strhtune)
    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'STRHOPT',NF_INT,1,strhopt)

    !
    ! Flags
    !
    istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'UBCFLG', NF_INT,1,ubcflg)
    istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'VBCFLG', NF_INT,1,vbcflg)
    istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'WBCFLG', NF_INT,1,wbcflg)
    istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'PTBCFLG',NF_INT,1,ptbcflg)
    istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'PRBCFLG',NF_INT,1,prbcflg)
    istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QVBCFLG',NF_INT,1,qvbcflg)

    istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'nscalar',NF_INT,1,nscalarout)

    IF (P_QC > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QCBCFLG',NF_INT,1,qscalarbcflg(P_QC))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QCBCFLG',NF_INT,1,0)
    END IF

    IF (P_QR > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QRBCFLG',NF_INT,1,qscalarbcflg(P_QR))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QRBCFLG',NF_INT,1,0)
    END IF

    IF (P_QI > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QIBCFLG',NF_INT,1,qscalarbcflg(P_QI))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QIBCFLG',NF_INT,1,0)
    END IF

    IF (P_QS > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QSBCFLG',NF_INT,1,qscalarbcflg(P_QS))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QSBCFLG',NF_INT,1,0)
    END IF

    IF (P_QH > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QHBCFLG',NF_INT,1,qscalarbcflg(P_QH))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QHBCFLG',NF_INT,1,0)
    END IF

    IF (P_QG > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QGBCFLG',NF_INT,1,qscalarbcflg(P_QG))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'QGBCFLG',NF_INT,1,0)
    END IF

    IF (P_NC > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NCBCFLG',NF_INT,1,qscalarbcflg(P_NC))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NCBCFLG',NF_INT,1,0)
    END IF

    IF (P_NR > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NRBCFLG',NF_INT,1,qscalarbcflg(P_NR))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NRBCFLG',NF_INT,1,0)
    END IF

    IF (P_NI > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NIBCFLG',NF_INT,1,qscalarbcflg(P_NI))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NIBCFLG',NF_INT,1,0)
    END IF

    IF (P_NS > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NSBCFLG',NF_INT,1,qscalarbcflg(P_NS))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NSBCFLG',NF_INT,1,0)
    END IF

    IF (P_NG > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NGBCFLG',NF_INT,1,qscalarbcflg(P_NG))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NGBCFLG',NF_INT,1,0)
    END IF

    IF (P_NH > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NHBCFLG',NF_INT,1,qscalarbcflg(P_NH))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NHBCFLG',NF_INT,1,0)
    END IF

    IF (P_ZR > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZRBCFLG',NF_INT,1,qscalarbcflg(P_ZR))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZRBCFLG',NF_INT,1,0)
    END IF

    IF (P_ZI > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZIBCFLG',NF_INT,1,qscalarbcflg(P_ZI))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZIBCFLG',NF_INT,1,0)
    END IF

    IF (P_ZS > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZSBCFLG',NF_INT,1,qscalarbcflg(P_ZS))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZSBCFLG',NF_INT,1,0)
    END IF

    IF (P_ZG > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZGBCFLG',NF_INT,1,qscalarbcflg(P_ZG))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZGBCFLG',NF_INT,1,0)
    END IF

    IF (P_ZH > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZHBCFLG',NF_INT,1,qscalarbcflg(P_ZH))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'ZHBCFLG',NF_INT,1,0)
    END IF

    IF (P_CC > 0) THEN
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'CCBCFLG',NF_INT,1,qscalarbcflg(P_CC))
    ELSE
      istatus = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'CCBCFLG',NF_INT,1,0)
    END IF

    ! do not fill, will set values explicitly later. Improve performance

    istatus = NF_SET_FILL(ncid,NF_NOFILL,oldfillmode)

!-----------------------------------------------------------------------
!
! Define variable arrays
!
!-----------------------------------------------------------------------

    istatus = NF_DEF_VAR(ncid,'CTIME',NF_CHAR,2,(/dims_id,dimt_id/),varid)
    CALL net_define_var_meta(ncid,varid,'REAL',bdymeta%ctime)

    IF (ubcflg /= 0) THEN
      istatus = NF_DEF_VAR(ncid,'U',NF_FLOAT,4,                         &
                         (/dimx_id,dimys_id,dimzs_id,dimt_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',bdymeta%u)
    END IF

    IF (vbcflg /= 0) THEN
      istatus = NF_DEF_VAR(ncid,'V',NF_FLOAT,4,                         &
                         (/dimxs_id,dimy_id,dimzs_id,dimt_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',bdymeta%v)
    END IF

    IF (wbcflg /= 0) THEN
      istatus = NF_DEF_VAR(ncid,'W',NF_FLOAT,4,                      &
                         (/dimxs_id,dimys_id,dimz_id,dimt_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',bdymeta%w)
    END IF

    IF (ptbcflg /= 0) THEN
      istatus = NF_DEF_VAR(ncid,'PT',NF_FLOAT,4,                      &
                         (/dimxs_id,dimys_id,dimzs_id,dimt_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',bdymeta%pt)
    END IF

    IF (prbcflg /= 0) THEN
      istatus = NF_DEF_VAR(ncid,'P',NF_FLOAT,4,                      &
                         (/dimxs_id,dimys_id,dimzs_id,dimt_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',bdymeta%p)
    END IF

    IF (qvbcflg /= 0) THEN
      istatus = NF_DEF_VAR(ncid,'QV',NF_FLOAT,4,                        &
                         (/dimxs_id,dimys_id,dimzs_id,dimt_id/),varid)
      CALL net_define_var_meta(ncid,varid,'REAL',bdymeta%qv)
    END IF

    DO nq = 1,nscalar
      IF (qscalarbcflg(nq) > 0) THEN
        istatus = NF_DEF_VAR(ncid,upcase(TRIM(qnames(nq))),NF_FLOAT,4,  &
                         (/dimxs_id,dimys_id,dimzs_id,dimt_id/),varid)
        CALL net_define_var_meta(ncid,varid,'REAL',bdymeta%qscalar(nq))
      END IF
    END DO

!-----------------------------------------------------------------------
!
! End NetCDF file DEFINE mode
!
!-----------------------------------------------------------------------

    istatus = NF_ENDDEF(ncid)
    CALL net_check_error(istatus,'net_define_exbc')

  END IF

  istatus = NF_INQ_VARID(ncid,'CTIME',varid)
  CALL net_check_error(istatus,'net_define_exbc')

  istatus = NF_PUT_VARA_TEXT(ncid,varid,(/1,itime/),(/ctime_len,1/),ctime)
  CALL net_check_error(istatus,'net_define_exbc')

  RETURN
END SUBROUTINE net_define_exbc
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_get_exbc                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_get_exbc(ncid,nx,ny,nz,itime,dx,dy,dz,                   &
             dzmin,strhopt,zrefsfc,dlayer1,dlayer2,zflat,strhtune,      &
             mapproj,sclfct,trulat1,trulat2,trulon,ctrlat,ctrlon,       &
             ubcflg,vbcflg,wbcflg,ptbcflg,prbcflg,qvbcflg,              &
             nscalarin,nqscalarbcflg,ctime,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Extract NetCDF file dimensions and attributes from ARPS boundary
!     data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (08/20/2004)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER, INTENT(IN)  :: ncid, itime
  INTEGER, INTENT(OUT) :: nx, ny, nz
  REAL,    INTENT(OUT) :: dx, dy, dz
  REAL,    INTENT(OUT) :: dzmin, strhtune
  INTEGER, INTENT(OUT) :: strhopt
  REAL,    INTENT(OUT) :: zrefsfc, dlayer1, dlayer2, zflat
  INTEGER, INTENT(OUT) :: mapproj
  REAL,    INTENT(OUT) :: sclfct
  REAL,    INTENT(OUT) :: trulat1, trulat2, trulon
  REAL,    INTENT(OUT) :: ctrlat,  ctrlon
  INTEGER, INTENT(OUT) :: ubcflg,   vbcflg,   wbcflg
  INTEGER, INTENT(OUT) :: ptbcflg,  prbcflg, qvbcflg
  INTEGER, INTENT(OUT) :: nscalarin
  INTEGER, INTENT(OUT) :: nqscalarbcflg(nscalar)

  CHARACTER(LEN=15), INTENT(OUT) :: ctime

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dimx_id,dimy_id,dimz_id,dimt_id,dims_id
  INTEGER :: lenstr

  INTEGER :: nq

  CHARACTER(LEN=4) :: upcase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Get dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_INQ_DIMID(ncid,'x_stag',dimx_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_exbc')
  istatus = NF_INQ_DIMLEN(ncid,dimx_id,nx)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_exbc')

  istatus = NF_INQ_DIMID(ncid,'y_stag',dimy_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_exbc')
  istatus = NF_INQ_DIMLEN(ncid,dimy_id,ny)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_exbc')

  istatus = NF_INQ_DIMID(ncid,'z_stag',dimz_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_exbc')
  istatus = NF_INQ_DIMLEN(ncid,dimz_id,nz)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_exbc')

!  istatus = NF_INQ_DIMID(ncid,'Time',dimt_id)
!  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_exbc')
!  istatus = NF_INQ_DIMLEN(ncid,dimt_id,nt)
!  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_exbc')

  istatus = NF_INQ_DIMID(ncid,'CtimeStrLen',dims_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_exbc')
  istatus = NF_INQ_DIMLEN(ncid,dims_id,lenstr)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_exbc')

!-----------------------------------------------------------------------
!
! Get global attributes
!
!-----------------------------------------------------------------------

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DX',dx)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:dx in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DY',dy)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:dy in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DZ',dz)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:dz in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DZMIN',dzmin)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:dzmin in net_get_exbc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'STRHOPT',strhopt)
  CALL net_check_error(istatus,'NF_GET_ATT_INT:strhopt in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'STRHTUNE',strhtune)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:strhtune in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'ZREFSFC',zrefsfc)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:zresfc in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DLAYER1',dlayer1)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:dlayer1 in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DLAYER2',dlayer2)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:dlayer2 in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'ZFLAT',zflat)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:zflat in net_get_exbc')

  !
  ! Map projection
  !
  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'MAPPROJ',mapproj)
  CALL net_check_error(istatus,'NF_GET_ATT_INT:mapproj in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'SCLFCT',sclfct)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:sclfct in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT1',trulat1)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:trulat1 in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT2',trulat2)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:trulat2 in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELON',trulon)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:trulon in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLAT',ctrlat)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:ctrlat in net_get_exbc')

  istatus = NF_GET_ATT_REAL(ncid,NF_GLOBAL,'CTRLON',ctrlon)
  CALL net_check_error(istatus,'NF_GET_ATT_REAL:ctrlon in net_get_exbc')

  !
  ! Flags
  !
  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'UBCFLG',ubcflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT:ubcflg in net_get_exbc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'VBCFLG',vbcflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT:vbcflg in net_get_exbc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'WBCFLG',wbcflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT:wbcflg in net_get_exbc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'PTBCFLG',ptbcflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT:ptbcflg in net_get_exbc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'PRBCFLG',prbcflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT:prbcflg in net_get_exbc')

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'QVBCFLG',qvbcflg)
  CALL net_check_error(istatus,'NF_GET_ATT_INT:qvbcflg in net_get_exbc')

  nscalarin = 0
  nqscalarbcflg(:) = 0

  istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'nscalar',nscalarin)
  IF (istatus == NF_NOERR) THEN   ! new version

    DO nq = 1,nscalar
      istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,                          &
                    TRIM(upcase(qnames(nq)))//'BCFLG',nqscalarbcflg(nq))
      CALL net_check_error(istatus,'NF_GET_ATT_INT:qgbcflg in net_get_exbc')
    END DO

  ELSE                            ! Maybe old version

    IF (P_QC > 0) THEN
      istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'QCBCFLG',nqscalarbcflg(P_QC))
      CALL net_check_error(istatus,'NF_GET_ATT_INT:qcbcflg in net_get_exbc')
      nscalarin = nscalarin + 1
    END IF

    IF (P_QR > 0) THEN
      istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'QRBCFLG',nqscalarbcflg(P_QR))
      CALL net_check_error(istatus,'NF_GET_ATT_INT:qcbcflg in net_get_exbc')
      nscalarin = nscalarin + 1
    END IF

    IF (P_QI > 0) THEN
      istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'QIBCFLG',nqscalarbcflg(P_QI))
      CALL net_check_error(istatus,'NF_GET_ATT_INT:qcbcflg in net_get_exbc')
      nscalarin = nscalarin + 1
    END IF

    IF (P_QS > 0) THEN
      istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'QSBCFLG',nqscalarbcflg(P_QS))
      CALL net_check_error(istatus,'NF_GET_ATT_INT:qcbcflg in net_get_exbc')
      nscalarin = nscalarin + 1
    END IF

    IF (P_QH > 0) THEN
      istatus = NF_GET_ATT_INT(ncid,NF_GLOBAL,'QHBCFLG',nqscalarbcflg(P_QH))
      CALL net_check_error(istatus,'NF_GET_ATT_INT:qcbcflg in net_get_exbc')
      nscalarin = nscalarin + 1
    END IF

  END IF


  istatus = NF_INQ_VARID(ncid,'CTIME',varid)
  CALL net_check_error(istatus,'NF_INQ_VARID:ctime in net_get_exbc')
  istatus = NF_GET_VARA_TEXT(ncid,varid,(/1,itime/),(/lenstr,1/),ctime)
  CALL net_check_error(istatus,'NF_GET_VARA_TEXT:ctime in net_get_exbc')

  RETURN
END SUBROUTINE net_get_exbc
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_define_onevar              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_define_onevar(ncid,nx,ny,nz,varname,varlongname,varunits, &
                             istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Define one variable in NetCDF file.
!     The netCDF file should be in DATA mode after this call.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (06/13/2005)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,      INTENT(IN) :: ncid
  INTEGER,      INTENT(IN) :: nx, ny, nz
  CHARACTER(*), INTENT(IN) :: varname, varlongname, varunits

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER           :: varid
  INTEGER           :: dimx_id, dimy_id, dimz_id
  INTEGER           :: lenstr

  INTEGER           :: oldfillmode

  CHARACTER(LEN=80) :: tmpstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Define dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_DEF_DIM(ncid,'x',nx,dimx_id)
  CALL net_check_error(istatus,'net_define_onevar')

  istatus = NF_DEF_DIM(ncid,'y',ny,dimy_id)
  CALL net_check_error(istatus,'net_define_onevar')

  istatus = NF_DEF_DIM(ncid,'z',nz,dimz_id)
  CALL net_check_error(istatus,'net_define_onevar')

!-----------------------------------------------------------------------
!
! Define global attributes
!
!-----------------------------------------------------------------------

  ! do not fill, will set values explicitly later. Improve performance

  istatus = NF_SET_FILL(ncid,NF_NOFILL,oldfillmode)

!-----------------------------------------------------------------------
!
! Define variable arrays
!
!-----------------------------------------------------------------------

  istatus = NF_DEF_VAR(ncid,varname,NF_FLOAT,3,                      &
                         (/dimx_id,dimy_id,dimz_id/),varid)

  lenstr = LEN_TRIM(varname)
  istatus = NF_PUT_ATT_TEXT(ncid,varid,'standard_name',lenstr,varname)

  lenstr = LEN_TRIM(varlongname)
  istatus = NF_PUT_ATT_TEXT(ncid,varid,'long_name',lenstr,varlongname)

  lenstr = LEN_TRIM(varunits)
  istatus = NF_PUT_ATT_TEXT(ncid,varid,'units',lenstr,varunits)


!-----------------------------------------------------------------------
!
! End NetCDF file DEFINE mode
!
!-----------------------------------------------------------------------

  istatus = NF_ENDDEF(ncid)
  CALL net_check_error(istatus,'net_define_onevar')

  RETURN
END SUBROUTINE net_define_onevar
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_get_onevar                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_get_onevar(ncid,nx,ny,nz,varname,varlongname,varunits,   &
                          istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!     Extract NetCDF file dimensions and variable attributes
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang (06/13/2005)
!
!  MODIFIED HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: ncid
  CHARACTER(*), INTENT(IN)  :: varname
  INTEGER,      INTENT(OUT) :: nx, ny, nz
  CHARACTER(*), INTENT(OUT) :: varlongname, varunits

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Included files
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dimx_id,dimy_id,dimz_id

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Get dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_INQ_DIMID(ncid,'x',dimx_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_onevar')
  istatus = NF_INQ_DIMLEN(ncid,dimx_id,nx)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_onevar')

  istatus = NF_INQ_DIMID(ncid,'y',dimy_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_onevar')
  istatus = NF_INQ_DIMLEN(ncid,dimy_id,ny)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_onevar')

  istatus = NF_INQ_DIMID(ncid,'z',dimz_id)
  CALL net_check_error(istatus,'NF_INQ_DIMID in net_get_onevar')
  istatus = NF_INQ_DIMLEN(ncid,dimz_id,nz)
  CALL net_check_error(istatus,'NF_INQ_DIMLEN in net_get_onevar')

!-----------------------------------------------------------------------
!
! Get variable attributes
!
!-----------------------------------------------------------------------

  istatus = NF_INQ_VARID(ncid,varname,varid)
  IF (istatus == NF_ENOTVAR) THEN  ! variable not found
    WRITE(6,'(1x,3a)') 'Variable <',trim(varname),'> is NOT FOUND.'
    istatus = 9999
    RETURN
  END IF
  !CALL net_check_error(istatus,'NF_INQ_VARID in net_get_onevar')

  istatus = NF_GET_ATT_TEXT(ncid,varid,'long_name',varlongname)
  CALL net_check_error(istatus,'NF_GET_ATT_TEXT in net_get_onevar')

  istatus = NF_GET_ATT_TEXT(ncid,varid,'units',varunits)
  CALL net_check_error(istatus,'NF_GET_ATT_TEXT in net_get_onevar')

  RETURN
END SUBROUTINE net_get_onevar
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_get_unlimit_size           ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE net_get_unlimit_size( filename, no_times, istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Get the size of unlimitted dimension in the file
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (04/20/2010)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(OUT) :: no_times
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  INTEGER :: ncid
  INTEGER :: dimid
  INTEGER :: timelen

  INCLUDE 'netcdf.inc'

  LOGICAL :: fexists

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus   = 0
  no_times  = 0

  INQUIRE(FILE = TRIM(filename), EXIST = fexists)

  IF (fexists) THEN
    istatus = NF_OPEN(TRIM(filename),NF_NOWRITE,ncid)
    CALL net_check_error(istatus,'netopen')

    istatus = NF_INQ_DIMID(ncid,'Time',dimid)

    IF (istatus == NF_NOERR) THEN
      istatus = NF_INQ_DIMLEN(ncid,dimid,timelen)
      CALL net_check_error(istatus,'NF_INQ_DIMLEN in get_ncd_frames_per_outfile')

      IF( timelen >= 1) THEN
        no_times = timelen
      ELSE
        WRITE(6,'(1x,3a,I2)') 'ERROR: The unlimited dimension in the file ',  &
                              TRIM(filename),' is bad, timelen = ',timelen, '.'
        istatus = -3
      END IF
    ELSE
      istatus = -2
    END IF

    istatus = NF_CLOSE( ncid )

  ELSE
    WRITE(6,'(2a)') 'File not found: ', filename
    istatus = -1
  END IF

  RETURN
END SUBROUTINE net_get_unlimit_size
