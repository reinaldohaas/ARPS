!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE calculate_hx             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE calculate_hx(nx,ny,nz,nzsoil,nstyps,timsnd,tmstrln,          &
            u_ct,v_ct,w_ct,wcont_ct,ptprt_ct,pprt_ct,                   &
            qv_ct,qscalar_ct,tke_ct,                                    &
            ubar_ct,vbar_ct,wbar_ct,ptbar_ct,pbar_ct,                   &
            rhobar_ct,qvbar_ct,kmh_ct,kmv_ct,                           &
            x_ct,y_ct,z_ct,zp_ct,zpsoil_ct,hterain_ct, mapfct_ct,       &
            soiltyp_ct,stypfrct_ct,vegtyp_ct,lai_ct,roufns_ct,veg_ct,   &
            tsoil_ct,qsoil_ct,wetcanp_ct,snowdpth_ct,qvsfc_ct,          &
            raing_ct,rainc_ct,prcrate_ct,                               &
            radfrc_ct,radsw_ct,rnflx_ct, radswnet_ct,radlwin_ct,        &
            usflx_ct,vsflx_ct,ptsflx_ct,qvsflx_ct,rad_file_state,       &
            sfc_file_state,snd_file_state,pro_file_state                &
            )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  calculate h(x) (radar data, surface data and other data)
!
!-----------------------------------------------------------------------
!  AUTHOR: Shizhang Wang
!    09/02/2010
!
!  This the program is designed according to SUBROUTINE enkfda, which
!  is designed and modified by:
!
!  AUTHOR: Mingjing Tong
!    09/05/2006.
!
!  MODIFIED HISTORY:
!    07/04/2007     Youngsun Jung
!       Add arps array module and polarimetric radar data
!    02/07/2009     Ting Lei
!       Added sfc enkf part
!    11/14/2009     Youngsun Jung
!       Added user defined radar
!-----------------------------------------------------------------------
!
!  INPUT/OUTPUT:
!
!    nx,ny,nz The dimension of data arrays
!-----------------------------------------------------------------------
  USE ARPSARRAY
  USE ARPS3DARRAY
  USE global_paraest
  USE define_observation_conventional
  USE PRIMARRAY
  USE module_mpi_enkf

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'arpsenkf.inc'
  INCLUDE 'mp.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz, tmstrln
  INTEGER :: nzsoil
  INTEGER :: nstyps
  CHARACTER (LEN=10) :: timsnd
  INTEGER :: istatus
  INTEGER :: i,j,k,n
  LOGICAL :: rad_file_state,sfc_file_state,snd_file_state,pro_file_state
!
!-----------------------------------------------------------------------
!
!  Input Variables.
!
!-----------------------------------------------------------------------
!
  REAL    :: x_ct(nx)
  REAL    :: y_ct(ny)
  REAL    :: z_ct(nz)
  REAL    :: mapfct_ct(nx,ny,8)
  REAL    :: zp_ct(nx,ny,nz)
  REAL    :: zpsoil_ct(nx,ny,nzsoil)
  REAL    :: hterain_ct(nx,ny)
  REAL    :: u_ct(nx,ny,nz)
  REAL    :: v_ct(nx,ny,nz)
  REAL    :: w_ct(nx,ny,nz)
  REAL    :: qv_ct(nx,ny,nz)
  REAL    :: wcont_ct(nx,ny,nz)
  REAL    :: ptprt_ct(nx,ny,nz)
  REAL    :: pprt_ct(nx,ny,nz)
  REAL    :: qscalar_ct(nx,ny,nz,1,nscalar)
  REAL    :: kmh_ct(nx,ny,nz)
  REAL    :: kmv_ct(nx,ny,nz)
  REAL    :: ubar_ct(nx,ny,nz)
  REAL    :: vbar_ct(nx,ny,nz)
  REAL    :: wbar_ct(nx,ny,nz)
  REAL    :: ptbar_ct(nx,ny,nz)
  REAL    :: pbar_ct(nx,ny,nz)
  REAL    :: rhobar_ct(nx,ny,nz)
  REAL    :: qvbar_ct(nx,ny,nz)
  INTEGER :: soiltyp_ct(nx,ny,nstyps)
  REAL    :: stypfrct_ct(nx,ny,nstyps)
  INTEGER :: vegtyp_ct(nx,ny)
  REAL    :: lai_ct(nx,ny)
  REAL    :: roufns_ct(nx,ny)
  REAL    :: veg_ct(nx,ny)
  REAL    :: tke_ct(nx,ny,nz)
  REAL    :: tsoil_ct(nx,ny,nzsoil,0:nstyps)
  REAL    :: qsoil_ct(nx,ny,nzsoil,0:nstyp)
  REAL    :: wetcanp_ct(nx,ny,0:nstyps)
  REAL    :: snowdpth_ct(nx,ny)
  REAL    :: qvsfc_ct(nx,ny,0:nstyps)
  REAL    :: raing_ct(nx,ny)
  REAL    :: rainc_ct(nx,ny)
  REAL    :: prcrate_ct(nx,ny,4)
  REAL    :: radfrc_ct(nx,ny,nz)
  REAL    :: radsw_ct(nx,ny)
  REAL    :: rnflx_ct(nx,ny)
  REAL    :: radswnet_ct(nx,ny)
  REAL    :: radlwin_ct(nx,ny)
  REAL    :: usflx_ct(nx,ny)
  REAL    :: vsflx_ct(nx,ny)
  REAL    :: ptsflx_ct(nx,ny)
  REAL    :: qvsflx_ct(nx,ny)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(rfopt == 0 .OR. rfopt == 10) CALL reflktb

  CALL allocateARPSarray(nx,ny,nz,nzsoil,nstyps)   ! Grid, base and static arrays
  CALL allocateARPS3Darray(nx,ny,nz,nzsoil,nstyps,nscalar)
  CALL allocateGridArray(nx,ny,nz)                 ! grid on mass grid points
  CALL allocateTruthArray(nx,ny,nz,nscalar)        ! Control run? 6+nscalar
  CALL allocateObsArray(nx,ny,nz,numvar)

  IF(paraestopt > 0) THEN
     CALL allocate_paraestArray(nen,paranum)          ! parameter estimates?
  ENDIF

  CALL allocateBarray(nx,ny,nz)                    ! bu,v,w,s, weightu,v,w,s

!
!-----------------------------------------------------------------------
!
!  Initial the value of arrays in ARPSarray
!
!-----------------------------------------------------------------------
!
  x(:)=x_ct(:)
  y(:)=y_ct(:)
  z(:)=z_ct(:)
  mapfct(:,:,:)=mapfct_ct(:,:,:)
  zp(:,:,:)=zp_ct(:,:,:)
  zpsoil(:,:,:)=zpsoil_ct(:,:,:)
  !uprt(:,:,:)=u_ct-ubar_ct(:,:,:)
  !vprt(:,:,:)=v_ct-vbar_ct(:,:,:)
  !wprt(:,:,:)=w_ct-wbar_ct(:,:,:)
  !qvprt(:,:,:)=qv_ct-qvbar_ct(:,:,:)
  kmh(:,:,:)=kmh_ct(:,:,:)
  kmv(:,:,:)=kmv_ct(:,:,:)
  ubar(:,:,:)=ubar_ct(:,:,:)
  vbar(:,:,:)=vbar_ct(:,:,:)
  wbar(:,:,:)=wbar_ct(:,:,:)
  ptbar(:,:,:)=ptbar_ct(:,:,:)
  pbar(:,:,:)=pbar_ct(:,:,:)
  rhobar(:,:,:)=rhobar_ct(:,:,:)
  qvbar(:,:,:)=qvbar_ct(:,:,:)
  soiltyp(:,:,:)=soiltyp_ct(:,:,:)
  stypfrct(:,:,:)=stypfrct_ct(:,:,:)
  vegtyp(:,:)=vegtyp_ct(:,:)
  lai(:,:)=lai_ct(:,:)
  roufns(:,:)=roufns_ct(:,:)
  veg(:,:)=veg_ct(:,:)

!
!-----------------------------------------------------------------------
!
!  Initial the value of arrays in ARPS4Darray
!
!-----------------------------------------------------------------------
!
  DO n=1,nen
  	u(:,:,:)=u_ct(:,:,:)
  	v(:,:,:)=v_ct(:,:,:)
  	w(:,:,:)=w_ct(:,:,:)
    ptprt(:,:,:)=ptprt_ct(:,:,:)
    pprt(:,:,:)=pprt_ct(:,:,:)
    qv(:,:,:)=qv_ct(:,:,:)
    DO i=1,nscalar
      qscalar(:,:,:,i)=qscalar_ct(:,:,:,1,i)
    ENDDO
    tke(:,:,:)=tke_ct(:,:,:)
    tsoil(:,:,:,:)=tsoil_ct(:,:,:,:)
    qsoil(:,:,:,:)=qsoil_ct(:,:,:,:)
    wetcanp(:,:,:)=wetcanp_ct(:,:,:)
    snowdpth(:,:)=snowdpth_ct(:,:)
    raing(:,:)=raing_ct(:,:)
    rainc(:,:)=rainc_ct(:,:)
    prcrate(:,:,:)=prcrate_ct(:,:,:)
    radfrc(:,:,:)=radfrc_ct(:,:,:)
    radsw(:,:)=radsw_ct(:,:)
    rnflx(:,:)=rnflx_ct(:,:)
    radswnet(:,:)=radswnet_ct(:,:)
    radlwin(:,:)=radlwin_ct(:,:)
    usflx(:,:)=usflx_ct(:,:)
    vsflx(:,:)=vsflx_ct(:,:)
    ptsflx(:,:)=ptsflx_ct(:,:)
    qvsflx(:,:)=qvsflx_ct(:,:)
  ENDDO
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  print*,'U'
!  print*,u
!  print*,'V'
!  print*,v
!  print*,'W'
!  print*,w
!  print*,'Zp'
!  print*,zp

  CALL enkfmpi_init_da( nx,ny, nz, nscalar, nen, r0hvr, r0hz, dx, dy,   &
                        lvldbg, istatus )

  CALL enkfmpi_state2ext1dx(nx,2,nx-2,xs,xs_ext,istatus)   ! They are time-independent
  CALL enkfmpi_state2ext1dx(nx,2,nx-1,x, x_ext, istatus)   ! do not need to be convert
  CALL enkfmpi_state2ext1dy(ny,2,ny-2,ys,ys_ext,istatus)   ! back.
  CALL enkfmpi_state2ext1dy(ny,2,ny-1,y, y_ext, istatus)
  CALL enkfmpi_state2ext3d(nx,ny,nz, zp,  zp_ext,  istatus)
  CALL enkfmpi_state2ext3d(nx,ny,nz, zps, zps_ext, istatus)
  CALL obsAtScalar(nx,ny,nz,u_ct,v_ct,w_ct,u,v,w)

!-----------------------------------------------------------------------
!
!  Output H(x) of conventional data
!
!-----------------------------------------------------------------------
!

  IF ((sfcassim == 1 .and. sfc_file_state) .OR.   &
      (sndassim == 1 .and. snd_file_state) .OR.   &
      (proassim == 1 .and. pro_file_state)        &
      ) THEN
    IF (myproc == 0) WRITE(6,'(1x,a)') '=== STEP H(x) output - conventional data'
     CALL calculate_hx_cnv(nx,ny,nz,timsnd,tmstrln,                     &
                           sfc_file_state,snd_file_state,pro_file_state)
!    CALL calculate_hx_cvn_sub(nx,ny,nz,timsnd,tmstrln,nzsoil,nstyps,sfc_file_state,snd_file_state,pro_file_state)
  ELSE
  	IF(.NOT. sfc_file_state) print*,'No sfc obs exists'
  	IF(sfcassim /= 1) print*,'sfc assimilation is disabled'
  	IF(.NOT. snd_file_state) print*,'No snd obs exists'
  	IF(sndassim /= 1) print*,'snd assimilation is disabled'
  	IF(.NOT. pro_file_state) print*,'No pro obs exists'
  	IF(proassim /= 1) print*,'pro assimilation is disabled'
  ENDIF

!
!-----------------------------------------------------------------------
!
!  Output H(x) of  radar data
!
!-----------------------------------------------------------------------
!
  IF(RADARDAOPT == 1 .and. rad_file_state  )THEN
    IF (myproc == 0) WRITE(6,'(1x,a)') '=== STEP H(x) output - radar data'
    CALL calculate_hx_rad(nx,ny,nz,timsnd,tmstrln)
  ELSE
  	IF(.NOT. rad_file_state ) print*,'No observation exists'
  	IF(RADARDAOPT /= 1 ) print*,'radar assimilation is disabled'
  ENDIF
!
!-----------------------------------------------------------------------
!
!  Deallocate work arrays
!
!-----------------------------------------------------------------------
!
  CALL deallocateARPSarray()   ! Grid, base and static arrays
  CALL deallocateARPS3Darray()
  CALL deallocateGridArray()                 ! grid on mass grid points
  CALL deallocateTruthArray()        ! Control run? 6+nscalar
  CALL deallocateObsArray()

  IF(paraestopt > 0) THEN
     CALL deallocate_paraestArray()          ! parameter estimates?
  ENDIF
  CALL deallocateBarray()                    ! bu,v,w,s, weightu,v,w,s

  !DEALLOCATE(rdrobs4D,obs_hgt,obs_rng,obs_elv)

  RETURN
END SUBROUTINE calculate_hx

!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE calculate_hx_rad            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
SUBROUTINE calculate_hx_rad(nx,ny,nz,timsnd,tmstrln)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Prior radar observations output
!
!-----------------------------------------------------------------------
!  AUTHOR: Shizhang Wang
!    09/02/2010
!
!  This the program is designed according to SUBROUTINE enkfradarda, which
!  is designed and modified by:
!
!  AUTHOR: Mingjing Tong
!    08/02/2006.
!
!  MODIFIED HISTORY:
!    07/04/2007     Youngsun Jung
!       Add arps array module and polarimetric radar data
!
!-----------------------------------------------------------------------

  USE ARPSARRAY
  USE ARPS3DARRAY
  USE radarobservations
  USE global_paraest
  USE radarInfo
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'arpsenkf.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz, tmstrln
  CHARACTER (LEN=10)  :: timsnd
  CHARACTER (LEN=256) :: rdrdtafile
  CHARACTER (LEN=3)   :: mem_num
  INTEGER :: temi01,temi02,temi03
  REAL    :: temr01,temr02,temr03
  INTEGER :: iradar, jradar,iflg
!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INTEGER :: i, j, k, idp
  INTEGER :: m, nkk, istatus, itilt
  INTEGER :: ni1,ni2,nj1,nj2,nk1,nk2
  INTEGER :: paraest
  REAL*8, ALLOCATABLE :: alpha(:,:,:,:)   ! Shape parameter

  CHARACTER(LEN=10) :: rdrnam
  INTEGER :: ntilt
  REAL    :: rdrlat,rdrlon,rdralt,radarx,radary, dazim,rngmin,rngmax
  INTEGER :: timeset, iyr, imon, iday, ihour, imin, isec

  INTEGER, ALLOCATABLE :: rdrdaopt(:)   ! See rdrdawopt/rdrdacopt in arpsenkf.input
  REAL,    ALLOCATABLE :: rdrthrshd(:)  ! See rdrthrshdw/rdrthrshdc in arpsenkf.input

  REAL, DIMENSION(:,:,:),     POINTER :: elvobs
  REAL, DIMENSION(:,:,:),     POINTER :: hgtoelv, rngoelv
  REAL, DIMENSION(:,:,:,:),   POINTER :: rdrobs
  REAL, DIMENSION(:,:,:,:,:), POINTER :: rdrobs_prior

  REAL, ALLOCATABLE :: rdrobs3d(:,:,:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  paraest= 0
  IF (SUM(paramest) > 0) THEN
    paraest = 1
  ENDIF
!-----------------------------------------------------------------
! Temporarily set the observation error to zero for hx calculation
!-----------------------------------------------------------------
!  rdrstd    Standard deviation of Gaussian error
!            (1)Vr  (2)Zhh  (3)Zvv  (4)Zdr  (5)Zdp  (6)Kdp
!            **<NOTE>** rdrstd(2)   STD in log domain (dbZ)

   rdrstd_vr = 0
   rdrstd_zhh= 0
   rdrstd_zvv= 0
   rdrstd_zdr= 0
   rdrstd_zdp= 0
   rdrstd_kdp= 0
   rdrstd_rhv= 0

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

  ALLOCATE(rdrdaopt(numvar), stat=istatus)
  CALL check_alloc_status(istatus, "enkfradarda:rdrdaopt")
  ALLOCATE(rdrthrshd(numvar), stat=istatus)
  CALL check_alloc_status(istatus, "enkfradarda:rdrthrshd")

  ALLOCATE(rdrobs3d(nx,ny,ntiltmax,numvar), stat=istatus)
  CALL check_alloc_status(istatus, "enkfradarda:rdrobs3d")

  ! First make sure parameter exptype is valid (radar data coordiates).
  IF (exptype /= 1 .AND. exptype /= 2) THEN
    WRITE(6,'(1x,a)') 'other type of radar data are not available now'
    WRITE(6,'(1x,a)') 'please set exptype <= 2'
    CALL arpsstop('ERROR: wrong exptype in enkradarda.',1)
  END IF

  IF(state2obsopt == 1)THEN
    iflg=1
  ELSE
    iflg=2
  ENDIF

!  IF(RADARDAOPT == 1)THEN

    IF(paraest > 0 .AND. corparaopt > 0)THEN
      IF(exptype == 1)THEN
        ALLOCATE(corflgrdr(nx,ny,nz,paranum,numvar), STAT=istatus)
                         ! correlation coefficient between parameter and Vr/Z
        CALL check_alloc_status(istatus, "enkfradarda:corflgrdr")
      ELSE IF(exptype == 2)THEN
        ALLOCATE(corflgrdr(nx,ny,ntiltmax,paranum,numvar), stat=istatus)
        CALL check_alloc_status(istatus, "enkfradarda:corflgrdr")
      END IF
    END IF

    DO k=1, nz-1
      DO j=1, ny-1
        DO i=1, nx-1
          tk(i,j,k)=(ptbar(i,j,k)+ptprt(i,j,k)) *                       &
                    (((pbar(i,j,k)+pprt(i,j,k))/p0)**rddcp)
        END DO
      END DO
    END DO

    SELECT CASE (rfopt)
    CASE (0)
      CALL reflec_ferrier_enkf(nx,ny,nz,rhobar,qscalar(:,:,:,P_QR),     &
                          qscalar(:,:,:,P_QS),                          &
                          qscalar(:,:,:,P_QH),                          &
                          tk,refdbz,rdr(:,:,:,2),zrs,zss,zhs)
    CASE (1:2)
      IF(numvar > 2) THEN
        CALL ZdrFromDualPol(nx,ny,nz,rhobar,qscalar,refdbz,             &
                            rdr(:,:,:,2),rdr(:,:,:,3),rdr(:,:,:,7),     &
                            rfopt,alpha)
        CALL Kdp(nx,ny,nz,rhobar,qscalar,rdr(:,:,:,6),rfopt,            &
                            wavelen,alpha)
      ELSE
        CALL ZhhFromDualPolCG(nx,ny,nz,rhobar,qscalar,                  &
                              refdbz,rdr(:,:,:,2),rfopt,wavelen,alpha)
      ENDIF
    CASE (10)
      CALL reflec_zhang(nx,ny,nz,rhobar,qscalar(:,:,:,P_QR),            &
                        qscalar(:,:,:,P_QS), qscalar(:,:,:,P_QH),       &
                        refdbz,rdr(:,:,:,2),zrs,zss,zhs)
    END SELECT

    DO m=1,nrdrused

      rdrdaopt(1)  = rdrdaopt1(m)
      rdrdaopt(2)  = rdrdaopt2(m)
      rdrthrshd(1) = rdrthd1(m)
      rdrthrshd(2) = rdrthd2(m)
      IF(numvar > 2) THEN
        rdrdaopt(3) = rdrdaopt3(m)
        rdrdaopt(4) = rdrdaopt4(m)
        rdrdaopt(5) = rdrdaopt5(m)
        rdrdaopt(6) = rdrdaopt6(m)
        rdrdaopt(7) = rdrdaopt7(m)
        rdrthrshd(3) = rdrthd3(m)
        rdrthrshd(4) = rdrthd4(m)
        rdrthrshd(5) = rdrthd5(m)
        rdrthrshd(6) = rdrthd6(m)
        rdrthrshd(7) = rdrthd7(m)
      ENDIF

      paraest=paramest(m)

      IF(exptype == 2)THEN
        !! Get radar parameteres  ! It is already called with obsdtaread, why call them again here - WYH?
        !SELECT CASE (ntwtype(m))
        !CASE (1:3)
        !  CALL init_radarPara(ntwtype(m),radarname(m),vcpmode(m),ntilt)
        !CASE (10)
        !  CALL assign_radarPara(U_nsweep,U_elvswp,U_ngate,U_gatesp,   &
        !                        U_beamwid,U_delaz,U_sradmul,          &
        !                        radarname(m),ntilt)
        !END SELECT

        radarn=0
        if(ntwtype(m) == 3) radarn= 1

      ENDIF

      elvobs  => obs_elv(:,:,:,m)
      hgtoelv => obs_hgt(:,:,:,m)
      rngoelv => obs_rng(:,:,:,m)

      !elvRDR  = obs_elv(:,m)
      !hgtoRDR = obs_hgt(:,:,:,m)
      !rngoRDR = obs_rng(:,:,:,m)

      !CALL getRDRpara(m)
      CALL getRDRpara(m,rdrnam,ntilt,                                   &
                      rdrlat,rdrlon,rdralt,radarx,radary,               &
                      dazim,rngmin,rngmax)

      IF(state2obsopt == 1)THEN
        CALL radarllxy(nx,ny,xs,ys,m,rdrlat,rdrlon,                     &
                       rdralt,radarx,radary,iradar,jradar)
      ENDIF

      CALL rdrgrdpara(nx,ny,nz,xs,ys,zps,m,                             &
                      azmsc,sfcrng,elvsc,rngsc,istatus)

      ni1 = 1
      ni2 = nx
      nj1 = 1
      nj2 = ny
      IF(exptype == 1)THEN
        nk1=1
        nk2=nz
        nkk=nz
      ELSE IF(exptype == 2)THEN
        nk1=1
        nk2=ntilt
        nkk=ntilt
      END IF

      print*,"calculate hx and output hx,obs and locations"

      IF(exptype == 1)THEN

        CALL state2obs_scgrd(iflg,nx,ny,nz,u,v,w,qscalar(:,:,:,P_QR),   &
                             qscalar(:,:,:,P_QS),qscalar(:,:,:,P_QH),   &
                             rhobar,azmsc,elvsc,rngsc,refdbz,rdr,       &
                             zrs,zss,zhs,rdrobs3D(:,:,:,:))
      ELSE IF(exptype == 2)THEN

        CALL state2obs_elv(iflg,nx,ny,nz,iradar,jradar,                 &
                           ntilt,rdrlat,rdrlon,rdralt,radarx,radary,    &
                           rngmin,rngmax,elvobs,hgtoelv,rngoelv,        &
                           u,v,w,qscalar,rdrobs3D)

      ELSE
        write(*,*)'additional data type are not available now'
        STOP
      ENDIF

      temi01=int(memid/100     )-int(memid/1000    )*10
      temi02=int(memid/10      )-int(memid/100     )*10
      temi03=int(memid/1       )-int(memid/10      )*10
      mem_num=char(48+temi01)//char(48+temi02)//char(48+temi03)

      CALL getRDRtime(m,timeset,iyr,imon,iday,ihour,imin,isec,istatus)

      !rdrnam = radarname(m)
      rdrdtafile=trim(rdrnam)//'_'//trim(hdmpfheader)//'_'//timsnd(1:tmstrln)//'_'//mem_num
      !print*,rdrdtafile
      !CALL dumprdrobs(nx,ny,nz,m,timsnd,tmstrln,curtim,rdrobs3D,rdrdtafile)


      CALL dumprdrobs(nx,ny,nz,nkk,numvar,rdrdtafile,rdrfmver,          &
                      hdmpfheader,exptype,0,vcpmode(m),ntwtype(m),      &
                      xs,ys,zps,                                        &
                      timeset,iyr,imon,iday,ihour,imin,isec,            &
                      rdrnam,rdrlat,rdrlon,radarx,radary,rdralt,        &
                      rngmin,rngmax,dazim,                              &
                      elvobs(2,2,:),hgtoelv,rngoelv,rdrobs3D,           &
                      timsnd,tmstrln,curtim,rdrobs3d,rdrobs3d,istatus)

      CALL deallocateRadarInfo()

    END DO ! DO m=1,nrdrused


    IF ( allocated(corflgrdr) ) DEALLOCATE (corflgrdr, STAT=istatus)

!  END IF ! IF(RADARDAOPT == 1)THEN

  DEALLOCATE (rdrdaopt, rdrthrshd)  ! Added by Y. Wang
  CALL deallocate_radar_da(istatus)
  DEALLOCATE (rdrobs3D)

  RETURN

END SUBROUTINE calculate_hx_rad


SUBROUTINE temp_change_vcpmode(curtim)

  USE radarInfo
  IMPLICIT NONE

  INCLUDE 'arpsenkf.inc'

  REAL     :: curtim
  INTEGER  :: iradar
  INTEGER  :: ielv

  ielv=mod(int(curtim/60),5)
  IF(ielv/=0) THEN
    PRINT*,"The ",ielv,"th elevation"
  ELSE
  	PRINT*,"The ",5,"th elevation"
  ENDIF
  DO iradar=1,nrdrused
    SELECT CASE(ielv)
      CASE(1)
        vcpmode(iradar)= 111
      CASE(2)
        vcpmode(iradar)= 112
      CASE(3)
        vcpmode(iradar)= 113
      CASE(4)
        vcpmode(iradar)= 114
      CASE(0)
        vcpmode(iradar)= 115
    END SELECT
  ENDDO

  print*,"covflgopt is temporary set to -1 for arps run"
  covflgopt = -1
  print*,"state2obsopt is temporary set to 2 for arps run"
  state2obsopt = 2
END SUBROUTINE temp_change_vcpmode

!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE calculate_hx_cnv            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
SUBROUTINE calculate_hx_cnv(nx,ny,nz,timsnd,tmstrln,                    &
                            sfc_file_state,snd_file_state,pro_file_state)

  USE ARPSARRAY
  USE ARPS3DARRAY
  USE define_observation_conventional

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'arpsenkf.inc'
  INCLUDE 'indtflg.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
!  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
! Conventional data
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: p(:,:,:)
  REAL, ALLOCATABLE :: theta(:,:,:)
  REAL, ALLOCATABLE :: td(:,:,:)
  REAL, ALLOCATABLE :: stnlat(:,:)
  REAL, ALLOCATABLE :: stnlon(:,:)
  REAL, ALLOCATABLE :: hterain(:,:)
  REAL, ALLOCATABLE :: wds(:,:,:)
  REAL, ALLOCATABLE :: wdd(:,:,:)

  REAL, ALLOCATABLE :: j1(:,:,:)
  REAL, ALLOCATABLE :: j2(:,:,:)
  REAL, ALLOCATABLE :: j3(:,:,:)
  REAL, ALLOCATABLE :: j3soil(:,:,:)
  REAL, ALLOCATABLE :: j3soilinv(:,:,:)

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz, ncyc, nkk
  CHARACTER (LEN=10) :: timsnd
  INTEGER :: nzsoil
  INTEGER :: nstyps
  INTEGER :: i,j,k,l,m,n,itilt,iflg

  INTEGER :: istatus           ! Staus indicator
  INTEGER :: istat

  INTEGER :: ireturn           ! Return status indicator
  INTEGER :: tmstrln
  INTEGER :: lfname
  INTEGER :: ni1, ni2, nj1, nj2, nk1, nk2
  INTEGER :: ntmp


  REAL, PARAMETER :: missvalue=-111.1
  REAL :: GASDEV
  REAL :: time
  REAL :: x3,y3,z3
  REAL :: rmsdvr, rmsdz, rmsdva

  REAL, ALLOCATABLE :: tem1dz(:), tem2dz(:), tem1(:,:,:)

  INTEGER :: CVN_ON
  LOGICAL :: sfc_file_state,snd_file_state,pro_file_state

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------
! Temporarily set the observation error to zero for hx calculation
!-----------------------------------------------------------------
!    err_sfc_uv   Standard observation error for u/v (m/s)
!    err_sfc_t    Standard observation error for t (C)
!    err_sfc_t_meso    Standard observation error for mesonet t (C)
!    err_sfc_td   Standard observation error for td (C)
!    err_sfc_p    Standard observation error for p (hPa)
!    err_sfc_pt   Standard observation error for pt (K)
!    err_sfc_qv   Standard observation error for qv (kg/kg)
!    err_sfc_rh   Standard observation error for rh

  err_sfc_uv = 0
  err_sfc_t = 0
  err_sfc_t_meso = 0
  err_sfc_td = 0
  err_sfc_p = 0
  err_sfc_pt = 0
  err_sfc_qv = 0
  err_sfc_rh = 0

!-----------------------------------------------------------------

  DO k=1, nz-1
    DO j=1, ny-1
      DO i=1, nx-1
        tk(i,j,k)=(ptbar(i,j,k)+ptprt(i,j,k)) *                         &
                  (((pbar(i,j,k)+pprt(i,j,k))/p0)**rddcp)
      END DO
    END DO
  END DO

!--------------------------------------------------------------------------
! Allocation arrays.
!--------------------------------------------------------------------------

  ALLOCATE(p(nx,ny,nz), stat=istatus)
  CALL check_alloc_status(istatus, "p")
  ALLOCATE(theta(nx,ny,nz), stat=istatus)
  CALL check_alloc_status(istatus, "theta")
  ALLOCATE(td(nx,ny,nz), stat=istatus)
  CALL check_alloc_status(istatus, "td")
  ALLOCATE(stnlat(nx,ny), stat=istatus)
  CALL check_alloc_status(istatus, "sfclat")
  ALLOCATE(hterain(nx,ny), stat=istatus)
  CALL check_alloc_status(istatus, "hterain")
  ALLOCATE(stnlon(nx,ny), stat=istatus)
  CALL check_alloc_status(istatus, "sfclon")
  ALLOCATE(wdd(nx,ny,nz), stat=istatus)
  CALL check_alloc_status(istatus, "wdd")
  ALLOCATE(wds(nx,ny,nz), stat=istatus)
  CALL check_alloc_status(istatus, "wds")

  ALLOCATE(j1(nx,ny,nz), stat=istatus)
  CALL check_alloc_status(istatus, "j1")
  ALLOCATE(j2(nx,ny,nz), stat=istatus)
  CALL check_alloc_status(istatus, "j2")
  ALLOCATE(j3(nx,ny,nz), stat=istatus)
  CALL check_alloc_status(istatus, "j3")
  ALLOCATE(j3soil(nx,ny,nzsoil), stat=istatus)
  CALL check_alloc_status(istatus, "j3soil")
  ALLOCATE(j3soilinv(nx,ny,nzsoil), stat=istatus)
  CALL check_alloc_status(istatus, "j3soilinv")

  ALLOCATE(tem1dz(nz),stat=istatus)
  ALLOCATE(tem2dz(nz),stat=istatus)
  ALLOCATE(tem1(nx,ny,nz),stat=istatus)

  CALL inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                          &
              hterain,mapfct,j1,j2,j3,j3soil,                           &
              j3soilinv,tem1dz,tem2dz,tem1)

  DEALLOCATE(tem1dz, tem2dz, tem1)

!--------------------------------------------------------------------------
! Get latitude and longitude of each grid points.
!--------------------------------------------------------------------------

  CALL xytoll(nx,ny,xs,ys,stnlat,stnlon)

!--------------------------------------------------------------------------
! Compute observed quantities.
!--------------------------------------------------------------------------
  DO k=1, nz-1
    DO j=1, ny-1
      DO i=1, nx-1
        p(i,j,k) = pbar(i,j,k) + pprt(i,j,k)
      END DO
    END DO
  END DO
  !
  ! Compute dew point temperature.
  !
  CALL getdew(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,p,tk,qv,td)
  !
  ! Dump out conventional data
  !
  CALL dump_conventional_observations(nx,ny,nz,zps,p,tk,td,u,v,         &
              stnlat,stnlon,xs,ys,hterain,                              &
              1,sfc_file_state,snd_file_state, pro_file_state)

  !print*,'pass 5'
  RETURN
END SUBROUTINE calculate_hx_cnv
