!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE CLOUD_CV                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE cloud_cv (nx,ny,nz,i4time,time,dirname,runname,              &
                     hdmpfmt,hdfcompr,lvldbg,                           &
                     xs,ys,zs,dx,dy,hterain,latgr,longr,                &
                     nxlg,nylg,xslg,yslg,                               &
                     p_3d,t_3d,rh_3d,                                   &
                     nobsng,indexsng,stnsng,isrcsng,csrcsng,xsng,ysng,  &
                     timesng,latsng,lonsng,hgtsng,                      &
                     kcloud,store_amt,store_hgt,                        &
                     ref_mos_3d,istat_radar,                            &
                     clouds_3d,cloud_ceiling,tem1,tem2d2,               &
                     istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform 3D fractional cloud cover analysis using SAO
!  cloud coverage observations, visible (VIS) and infrared (IR)
!  satellite measurements, and radar reflectivity data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  02/01/96
!
!  MODIFICATION HISTORY:
!
!  07/12/96  J. Zhang
!            Added remapping of the cloud cover fields onto the
!            ARPS grid.
!  03/14/97  J. Zhang
!            Cleaning up the code  and implemented for the official
!            arps4.2.4 version.
!  03/28/97  J. Zhang
!            Using the lifiting condensation level as the cloud
!            base when inserting radar reflectivity with no SAO
!            cloud base.
!  04/11/97  J. Zhang
!            Included dims.inc for nx,ny,nz.
!  04/15/97  J. Zhang
!            Change the size of the arrays istn, jstn, hstn, and
!            cstn from fixed (40) to adjustable (mx_sng)
!  08/06/97  J. Zhang
!            Change adascld24.inc to adascld25.inc
!  09/04/97  J. Zhang
!            Changed the radar echo thresholds for inserting clouds
!            from radar reflectivities.
!  09/10/97  J. Zhang
!            Calculate lifting condensation level in INICLDGRD,
!            and used by other routines (INSERTIR and INSERTRAD).
!            Change adascld25.inc to adascld26.inc
!  11/17/97  J. Zhang
!            Change the length of "runname" from fixed number
!            of letters (6) to flexible numbers.
!  11/18/97  J. Zhang
!            Added flag 'clddiag' in adascld26.inc for whether to
!            output the special cloud fields.
!  03/26/98  J. Zhang
!            Some modifications based on the official adas 4.3.3.
!  04/27/98  J. Zhang
!            Conform WRTCLDVAR to WRTVAR
!  05/05/98  J. Zhang
!            Abandoned the cloud grid, using the ARPS grid instead.
!  10/25/98  K. Brewster
!            Modified creation of file name for cloud observation
!            diagnostic file.
!  08/21/01  K. Brewster
!            Changed wrtvar calls to wrtvar1 to match new routine
!            for writing 3-D arrays to be read by arpsplt as an
!            arbitary array.
!  03/19/03  K. Brewster
!            Added code for new satellite calibration methods.
!            Streamlined code and changed some variable names to be
!            consistant with other ADAS code.
!
!  04/12/07  Y. Wang
!            Changed all wrtvar1 calls to wrtvar2 because wrtvar2
!            supports HDF format and netCDF fromat and it also can
!            join variable automatically in MPI mode.
!            Note that the dummy argument were added, hdmpfmt, hdfcompr.
!
!            Added dummy argument lvldbg for controlling of the
!            stardard outputs.
!
!-----------------------------------------------------------------------
!
!  INCLUDE: (from adas.inc)
!
!  mx_sng            ! max. possible # of surface obs.
!  max_cld_snd       ! max. possible # of stations
                     ! with cloud coverage reports.
!  refthr1           ! "significant" radar echo at lower levels
!  refthr2           ! "significant" radar echo at upper levels
!  hgtrefthr         ! height criteria for "significant" radar
                     ! echo thresholds
!  clddiag           ! flag for whether to output the special cloud
                     ! fields
!
!  INPUT:
!
!  nx,ny,nz          ! ARPS grid size
!
!  i4time            ! analysis time in seconds from 00:00 UTC
                     ! Jan. 1, 1960
!  time              ! analysis time in seconds from the model
                     ! initial time
!  dirname           ! directory for data dump
!
!  INPUT (for ARPS grid variables)
!   xs     (nx)      ! The x-coord. of the physical and
                     ! computational grid. Defined at p-point.
!   ys     (ny)      ! The y-coord. of the physical and
                     ! computational grid. Defined at p-point.
!   zs    (nx,ny,nz) ! The physical height coordinate defined at
                     ! p-point of the staggered grid.
!   temp_3d(nx,ny,nz)! the temperature field
!
!   hterain (nx,ny)     ! The height of the terrain (equivalent
!   dx, dy           ! grid spacing
!
!  INPUT for SAO obs.
!
!   stnsng (mx_sng)      ! station name of single-lvel data
!   isrcsng (mx_sng)     ! flag for stations used/not used
!   csrcsng (mx_sng)     ! names of sfc sources
!   timesng (mx_sng)     ! time for the observation.
!   latsng (mx_sng)      ! latitude of single-level data
!   lonsng (mx_sng)      ! longitude of single-level data
!   hgtsng (mx_sng)     ! height of single-level data
!   xsng (mx_sng)        ! x location of single-level data
!   ysng (mx_sng)        ! y location of single-level data
!
!   kcloud (mx_sng)      ! number of obs. cloud layers.
!   store_amt (mx_sng,5) ! cloud coverage (ea. layer,
                         ! ea. station).
!   store_hgt (mx_sng,5) ! height of obs. cloud layers.
!
!  LOCAL
!
!   cf_modelfg (nx,ny,nz) ! first guess cloud cover field
!   rh_modelfg (nx,ny,nz) ! first guess relative humidity field
!   t_sfc_k (nx,ny)       ! surface temperature field
!
!  OUTPUT:
!
!  istatus                ! The flag indicating process status.
!  clouds_3d (nx,ny,nz)   ! final gridded cloud cover analysis
!  cloud_ceiling (nx,ny)  ! cloud ceiling (M AGL)
!
!  LOCAL:
!  n_cld_snd         number of cloud soundings created
!  cld_snd (max_cld_snd,nz)   Cld snding obtained from SAO data.
!  wt_snd (max_cld_snd,nz)    weights assigned to cloud sounding
!                                obtained from SAO data.
!  cvr_max (nx,ny)           ! final column maximum cloud cover
!  cvr_tot (nx,ny)           ! final column total cloud cover
!  cloud_top (nx,ny)         ! cloud top (M ASL)
!
!  LOCAL for SAO:
!
!  ista_snd (max_cld_snd)    index of cloud sounding stations
!  i_snd (max_cld_snd)       i-locn of each cld snd stn in ADAS grid
!  j_snd (max_cld_snd)       j-locn of each cld snd stn in ADAS grid
!
!
!  LOCAL for gridded cloud cover analysis
!
!  cldcv (nx,ny,nz)     3D gridded fractional cloud cover analysis.
!                        (when input, it's the first guess field)
!  wtcldcv (nx,ny,nz)   wts given to gridded fracn cld cvr anx.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
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
  INCLUDE 'adas.inc'     ! ADAS parameters.
  INCLUDE 'adassat.inc'  ! ADAS satellite parameters.
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  INPUT: general
  INTEGER :: nx,ny,nz          ! ARPS grid size
  INTEGER :: nxlg,nylg         ! LARGE grid size (MPI)
  INTEGER :: i4time            ! analysis time in seconds from
                               ! 00:00 UTC Jan. 1, 1960
  REAL    :: time               ! analysis time in seconds from the
                               ! model initiali time
  CHARACTER (LEN=*) :: dirname     ! directory for data dump
  CHARACTER (LEN=*) :: runname     ! a string identify the run

  INTEGER, INTENT(IN) ::  hdmpfmt,hdfcompr,lvldbg
!
!  INPUT: model grid
!
  REAL :: dx,dy                ! ADAS grid spacing

  REAL :: xs    (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at p-point.
  REAL :: ys    (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at p-point.
  REAL :: xslg(nxlg)           ! Same as xs for large grid
  REAL :: yslg(nylg)           ! Same as ys for large grid
  REAL :: zs (nx,ny,nz)        ! The physical height coordinate defined at
                               ! p-point of the staggered grid.
  REAL :: t_3d(nx,ny,nz)       ! temperature (K)
  REAL :: p_3d(nx,ny,nz)
  REAL :: rh_3d(nx,ny,nz)
  REAL :: hterain (nx,ny)         ! The height of the terrain (equivalent
!
!  INPUT for SAO obs.
!
  INTEGER :: nobsng
  INTEGER:: indexsng(nobsng)    !  which processor "owns" the ob
  CHARACTER (LEN=5) :: stnsng(mx_sng)    ! station name of single-lvel data
  INTEGER :: isrcsng(mx_sng)    ! is station used?
  CHARACTER (LEN=8) :: csrcsng(mx_sng)! names of sfc sources
  INTEGER :: timesng(mx_sng)    ! time for the observation.
  REAL :: xsng(mx_sng)          ! x location of single-level data
  REAL :: ysng(mx_sng)          ! y location of single-level data
  REAL :: latsng(mx_sng)        ! latitude of single-level data
  REAL :: lonsng(mx_sng)        ! longitude of single-level data
  REAL :: hgtsng(mx_sng)       ! height of single-level data
!
  INTEGER :: kcloud(mx_sng)     ! number of cloud layers.
  CHARACTER (LEN=4) :: store_amt(mx_sng,5) ! cld coverage (ea.lyr, ea. stn).
  REAL :: store_hgt (mx_sng,5)  ! height of cloud layers.
!
!  INPUT: for radar
  INTEGER :: istat_radar        ! Valid radar reflectivity data?
  REAL :: ref_mos_3d(nx,ny,nz)  ! 3D radar reflectivity on model grid
                                ! may be modified by sate VIS data
!
!  OUTPUT:
  INTEGER :: istatus              ! Flag for process status
  REAL :: clouds_3d(nx,ny,nz)     ! Final cloud cover field
  REAL :: cloud_ceiling (nx,ny)   ! cloud ceiling (M AGL)
  REAL :: tem1(nx,ny,nz)
  REAL :: tem2d2(nx,ny,2)
!
!  LOCAL: for first guess fields
  REAL :: cf_modelfg (nx,ny,nz)   ! first guess cloud cover field
  REAL :: rh_modelfg (nx,ny,nz)   ! first guess relative humidity fld
  REAL :: z_lcl (nx,ny)           ! lifting condensation level(MSL)
!
!  LOCAL: for surface arrays
  REAL :: latgr (nx,ny)           ! latitude at each grid point
  REAL :: longr (nx,ny)           ! longitude at each grid point
  REAL :: pres_sfc_pa(nx,ny)      ! sfc pressure (Pascal)
  REAL :: t_gnd_k(nx,ny)          ! ground skin temp.
  REAL :: t_sfc_k (nx,ny)         ! surface temperature field
  REAL :: cvr_snow (nx,ny)        ! surface snow coverage
  REAL :: rland_frac (nx,ny)      ! land coverage fraction
!
!  LOCAL: for SAO cloud sounding
  INTEGER :: ista_snd (max_cld_snd)   ! seq. index of cld snd stn
  INTEGER :: i_snd (max_cld_snd)      ! i-loctn of cld snd
  INTEGER :: j_snd (max_cld_snd)      ! j-loctn of cld snd
  INTEGER :: n_cld_snd                ! # of cloud snds created
  REAL :: cld_snd (max_cld_snd,nz)    ! SAO cld snd
  REAL :: wt_snd (max_cld_snd,nz)     ! wgt for SAO cld snd
!
!  LOCAL for gridded cloud cover analysis
  REAL :: cldcv (nx,ny,nz)       ! 3D gridded frac cld cv analysis.
                               ! (it's also bkgrnd field when input)
  REAL :: wtcldcv (nx,ny,nz)     ! wts assgned to cld cvr analysis
                               ! (it's also bkgrnd field when input)
!  LOCAL: for solar positions
!
  REAL :: solar_alt(nx,ny)       ! solar altitude angle
  REAL :: solar_ha(nx,ny)        ! solar hour angle
  REAL :: solar_dec              ! solar declination angle
!
!  LOCAL: for satellite insertion
!
  INTEGER :: io_vis              ! Valid satellite VIS albedo data?
  INTEGER :: isatvis(nx,ny)      ! Index of sat for Mosaicked vis data
  REAL :: albedo(nx,ny)          ! vis. albedo derived from sat vis data
  REAL :: cloud_frac_vis_a(nx,ny)! cldcvr derived from vis albedo
  INTEGER :: istat_vis           ! status of inserting vis data

  INTEGER :: io_ir               ! Valid satellite bright. temp. data?
  INTEGER :: io_calib            ! Valid satellite calibration data?
  INTEGER :: isatir(nx,ny)       ! Index of sat for mosaicked IR sat data
  REAL :: irtemp(nx,ny,2)        ! Mosaicked 10.7-micron. temp.
                                 ! 1: temp 2: cold-filtered 10.7 temp
  REAL :: cldtop_m_irt(nx,ny)    ! Sat. cld top ht (m) from 10-micron data
  REAL :: cldtop_m (nx,ny)       ! cldtop height
  REAL :: cloud_top (nx,ny)      ! cloud top (M ASL)
  INTEGER :: istatus_ir          ! status of inserting IR data
!
  REAL :: sfc_sao_buffer         ! No clearing cloud by sat. below
                                 ! this ht.(m AGL) (hence letting
                                 ! SAOs dominate)
  PARAMETER (sfc_sao_buffer = 800.) ! meters, keep lower SAO cld from
                                    ! cleared out by sate. ir data
!
!  LOCAL: for radar data insertion
  REAL :: dbz_max_2d(nx,ny)       ! column max. radar refl. (dbZ)
  REAL :: cloud_base (nx,ny)      ! cloud base heights from SAO+sate
  REAL :: cloud_base_buf(nx,ny)   ! lowest cloud base heights within
                                  ! a searching radius
  LOGICAL :: l_unresolved(nx,ny)  ! no SAO+sate cld, yet radar echo is
                                  ! above the threshold
  REAL :: vis_rad_thdbz
  PARAMETER (vis_rad_thdbz=10.)! threshold define cloudy refl.
  REAL :: vis_rad_thcvr
  PARAMETER (vis_rad_thcvr=0.2)! VIS cld cvr threshold,
                               ! below this thresh. may cause
                               ! radar echo being cleared out.
!
!  LOCAL: intermidiate products
  REAL :: cvr_max (nx,ny)         ! column maximum cloud cover
  REAL :: cvr_tot (nx,ny)         ! final column total cloud cover
  REAL :: cvr_1d (nz)             ! cloud cover array for 1 grid col.
!
!  LOCAL: for post-process
  REAL :: thresh_cvr_base,thresh_cvr_top,thresh_cvr_ceiling
  PARAMETER (thresh_cvr_base = 0.1)
  PARAMETER (thresh_cvr_top  = 0.1)
  PARAMETER (thresh_cvr_ceiling = 0.65)

  REAL :: default_top,default_base,default_ceiling
  PARAMETER (default_base     = 0.0)
  PARAMETER (default_top      = 0.0)
  PARAMETER (default_ceiling  = 0.0)

!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k1,l,j1,j2,ibeg
  LOGICAL :: l_stn_clouds       ! Using SAO stns' cloud obs?
  PARAMETER (l_stn_clouds = .true.)
  LOGICAL :: cldprt
!
  REAL :: default_clear_cover
  PARAMETER (default_clear_cover=0.01)
!
  INTEGER :: istatus_sao,istatus_fg,istatus_barnes                      &
          ,istatus_rad,istatus_vis
  INTEGER :: i,j,k,imid,jmid
  INTEGER :: imidproc, jmidproc
  REAL :: grid_spacing_m
  REAL :: albmin,albmax,vcfmin,vcfmax

  CHARACTER (LEN=6) :: varid
  CHARACTER (LEN=20) :: varname
  CHARACTER (LEN=20) :: varunits

  INTEGER :: nstn,istn(mx_sng),jstn(mx_sng)
  REAL :: hstn(mx_sng)
  CHARACTER (LEN=5) :: cstn(mx_sng)

  INTEGER :: istat_col

  CHARACTER (LEN=6)   :: satnam
  CHARACTER (LEN=256) :: fname
  REAL :: latsat,lonsat
  INTEGER :: itime,isource,ifield,lrunnam,istatwrt
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,*)' Welcome to the ADAS Cloud Analysis'
  cldprt=(clddiag == 1)
  ibeg=MAX(1,(nx-10))
  istatus=0

  imid = nxlg/2     ! Added by Y. Wang for proper diagnostic outputs
  jmid = nylg/2
  imidproc = (imid-2) / (nx-3) + 1
  jmidproc = (jmid-2) / (ny-3) + 1

  IF (loc_x == imidproc) THEN
    imid = MOD((imid-2),(nx-3)) + 2   ! Local index for global middle point
  ELSE
    imid = -999
  END IF

  IF (loc_y == jmidproc) THEN
    jmid = MOD((jmid-2),(ny-3)) + 2   ! Local index for global middle point
  ELSE
    jmid = -999
  END IF

!
!-----------------------------------------------------------------------
!
!  Initialize cloud cover fields on the analysis grid
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) PRINT*,'Call background fields routine'
  CALL inicldgrd(nx,ny,nz,zs,                                           &
          bgqcopt,default_clear_cover,                                  &
          hterain,t_3d,p_3d,rh_3d,cf_modelfg,t_sfc_k,pres_sfc_pa,       &
          cldcv,wtcldcv,z_ref_lcl,z_lcl,rh_modelfg,                     &
          r_missing,tem1,istatus_fg)
!
  IF (cldprt .AND. loc_y == jmidproc .AND. loc_x == nproc_x) THEN
    WRITE(6,'(/a)') ' ==== CLOUD_CV cloud cover first guess===='
    WRITE(6,'(/1X,3X,11(1X,i4,1X))') (i,i=ibeg,nx)
    WRITE(6,'(1X,i3,11F6.1)')                                           &
         (k,(cf_modelfg(i,jmid,k),i=ibeg,nx),k=nz,1,-1)
  END IF
!
  IF (istatus_fg /= 1)THEN
    IF (myproc == 0)                                                    &
      WRITE(6,*)' Error in background field: Aborting cloud analysis'
    GO TO 999
  END IF
!
  IF (cloudopt == 2) THEN   ! Using radar data only
    IF (myproc == 0)                                                    &
      PRINT*,'## NOTE ##:  Using radar data only. LCL is used.'
    GO TO 330
  END IF
!
!-----------------------------------------------------------------------
!
!  Create cloud sounding from the SAO cloud observations.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,'(/a/)') ' Call Ingest/Insert SAO routines'

  n_cld_snd = 0
  IF (nobsng >= 1) THEN
!
    CALL insert_sao1 (nx,ny,nz,dx,dy,xs,ys,zs,hterain,t_sfc_k           &
              ,nxlg,nylg,xslg,yslg                                      &
!              ,cldcv,wtcldcv,t_3d,rh_modelfg                            &
              ,cldcv,wtcldcv,t_3d,cf_modelfg      & ! fixed by Guoqing Ge
              ,nobsng,indexsng,stnsng,isrcsng,csrcsng                   &
              ,xsng,ysng,ista_snd                                       &
              ,timesng,latsng,lonsng,hgtsng                             &
              ,kcloud,store_amt,store_hgt                               &
              ,l_stn_clouds,n_cld_snd,cld_snd,wt_snd                    &
              ,i_snd,j_snd                                              &
              ,istatus_sao)
!
    IF(istatus_sao /= 1) THEN
      IF (myproc == 0)                                                  &
        WRITE(6,*)' Error inserting SAO data: Aborting cloud analysis'
      GO TO 999
    END IF
  ELSE
    IF (myproc == 0)                                                    &
     WRITE(6,'(a,a,/a)')                                                &
     ' !!!WARNING!!! No surface cloud obs. available.',                 &
                     ' Skipping SAOs insertion'
  END IF
!
!-----------------------------------------------------------------------
!
!  Initialize clouds_3d array
!
!-----------------------------------------------------------------------
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        clouds_3d(i,j,k)=0.0
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Do grid analysis on SAO data using Barnes interpolation.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,'(/a)') ' Analyzing SAO cloud soundings'
!
  grid_spacing_m = dx

  CALL barnes_r5 (nx,ny,nz,cldcv,clouds_3d,wtcldcv,cf_modelfg           &
          ,l_stn_clouds,default_clear_cover,cld_snd,wt_snd              &
          ,grid_spacing_m,i_snd,j_snd,n_cld_snd                         &
          ,istatus_barnes)
!
  IF(cldprt .AND. loc_y == jmidproc .AND. loc_x == nproc_x) THEN
    WRITE(6,'(/a)') ' ==== CLOUD_CV cloud cover barnes===='
    WRITE(6,'(/1x,3x,11(1x,i4,1x))') (i,i=ibeg,nx)
    WRITE(6,'(1x,i3,11f6.1)') &
         (k,(clouds_3d(i,jmid,k),i=ibeg,nx),k=nz,1,-1)
  END IF

  IF (istatus_barnes /= 1) THEN
    IF (myproc == 0)                                                    &
      WRITE(6,*)' Error in Barnes interp, Aborting cld analysis...'
    GO TO 999
  END IF
!
!-----------------------------------------------------------------------
!
!  Find all stations within the model domain and print the cloud
!  soundings at these station points after the Barnes interpolation.
!
!-----------------------------------------------------------------------
!
  k=0
  DO i=1,n_cld_snd
    IF(i_snd(i) >= 1.AND.i_snd(i) <= nx                                 &
          .AND. j_snd(i) <= ny.AND.j_snd(i) >= 1) THEN
      k=k+1
      cstn(k)=stnsng(ista_snd(i))
      istn(k)=i_snd(i)
      jstn(k)=j_snd(i)
      hstn(k)=0.001*hgtsng(ista_snd(i))
    END IF
  END DO
  nstn=k

  IF (cld_files == 1 .AND. myproc == 0) THEN
    IF (mp_opt > 0) WRITE(6,*) 'CLD_FILES for processor 0 only'
    CALL gtlfnkey( runname, lrunnam )
    WRITE(fname,'(a,a)') runname(1:lrunnam),'.cld_brns_snd'
    OPEN (UNIT=14,FILE=trim(fname),STATUS='unknown')

    WRITE(14,'(1x,a,2i5/a,i3/)')                                       &
          ' Total/in-domain # of stations:',nobsng,nstn,               &
          ' Cloud observations used for deriving cloud soundings:',    &
          n_cld_snd

    i = nstn/24

    IF(i < 1) THEN

      WRITE(14,'(1x,a,1x,24a5)') 'hgt(m)',(cstn(k),k=1,nstn)
      WRITE(14,'(8x,24i5)') (istn(k),k=1,nstn)
      WRITE(14,'(8x,24i5)') (jstn(k),k=1,nstn)
      WRITE(14,'(8x,24f5.2)') (hstn(k),k=1,nstn)
      DO k1=1,nz
        WRITE(14,'(1x,F7.1,24F5.2)') zs(istn(1),jstn(1),k1)             &
                   ,(clouds_3d(istn(l),jstn(l),k1),l=1,nstn)
      END DO !k1
    ELSE

      DO k=1,i
        j1=(k-1)*24+1
        j2=k*24
        WRITE(14,'(1x,a,1x,a5)') 'hgt(m)',cstn(k)
        WRITE(14,'(8x,i5)') istn(k)
        WRITE(14,'(8x,i5)') jstn(k)
        WRITE(14,'(8x,f5.2)') hstn(k)

        DO k1=1,nz
          WRITE(14,'(1x,f7.1,24f5.2)') zs(istn(j1),jstn(j1),k1),        &
                    (clouds_3d(istn(l),jstn(l),k1),l=j1,j2)
        END DO !k1
      END DO !k

      j=i*24+1
      WRITE(14,'(1x,a,1x,24a5)') 'hgt(m)',(cstn(l),l=j,nstn)
      WRITE(14,'(8x,24i5)') (istn(l),l=j,nstn)
      WRITE(14,'(8x,24i5)') (jstn(l),l=j,nstn)
      WRITE(14,'(8x,24f5.2)') (hstn(l),l=j,nstn)
      DO k1=1,nz
        WRITE(14,'(1x,F7.1,24F5.2)') zs(istn(j),jstn(j),k1),            &
                    (clouds_3d(istn(l),jstn(l),k1),l=j,nstn)
      END DO !k1
    END IF
    CLOSE(14)
  END IF
!
  IF (cld_files == 1) THEN
    IF (myproc == 0) WRITE(6,*)  'Writing 3D cloud fraction after inserting SAOs.'
    varid='cldsao'
    varname='Clouds after SAO'
    varunits='Fraction'
    CALL wrtvar2(nx,ny,nz,clouds_3d,varid,varname,varunits,             &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in the surface snow coverage data.
!  Currently, the snow cover field is simply set to zero
!
!-----------------------------------------------------------------------
!
!  varname='snwcvr'
!  CALL READVAR(nx,ny,1, varname,time,cvr_snow, runname)

  DO j = 1,ny
    DO i = 1,nx
      cvr_snow(i,j) = 0.0
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Read in the land (water) coverage fraction data.
!  Currently, this is simply set to 1.
!
!-----------------------------------------------------------------------
!
!  varname='lndfrc'
!  CALL READVAR(nx,ny,1, varname,time,rland_frac, runname)

  DO j = 1,ny
    DO i = 1,nx
      rland_frac(i,j) = 1.0
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate solar altitude
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,*)  'Calculating solar position parameters.'
!
  DO j = 1,ny
    DO i = 1,nx
      CALL solar_pos( latgr(i,j),longr(i,j),                            &
                      i4time,solar_alt(i,j),solar_dec,solar_ha(i,j) )
    END DO
  END DO
!
! For the MPI case, the next output presents processor 0 instead of the
! entire domain.  This is harmless.
!
  IF (loc_x == imidproc .AND. loc_y == jmidproc)                        &
    WRITE(6,'(3x,2(a,I3),a,3(a,F12.5))') 'pt(',imid,',',jmid,'): ',     &
            'Solar alt = ',solar_alt(imid,jmid),                        &
            ', declin = ',solar_dec,', hour_angle = ',solar_ha(imid,jmid)
!
!-----------------------------------------------------------------------
!
!  READ IN SATELLITE VISIBLE DATA
!
!-----------------------------------------------------------------------
!
  ifield=1
  isource=1

  IF (myproc == 0) WRITE(6,*)   &
    'Getting satellite visible data remapped on the model grid'

  CALL vismosaic(nx,ny,mx_sat,nvisfiles,vis_fname,viscalname,           &
                     isatvis,albedo,tem2d2,                             &
                     satnamvis,latsatvis,lonsatvis,io_vis)

  IF (cld_files == 1) THEN
    varid='albedo'
    varname='Albedo'
    varunits='Fraction'
    CALL wrtvar2(nx,ny,1 ,albedo, varid,varname,varunits,               &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
  END IF

  IF(io_vis == 0) THEN
    IF (myproc == 0) WRITE(6,'(1x,a)') 'Read successful, calling process_vis'
    CALL process_vis (solar_alt,nx,ny,                                  &
                   cloud_frac_vis_a,albedo,r_missing,lvldbg,istat_vis)

    albmax = -999.0
    albmin =  999.0
    vcfmax = -999.0
    vcfmin =  999.0

    DO j=1,ny
      DO i=1,nx
        IF(albedo(i,j) >= 0.) THEN
            albmax=MAX(albmax,albedo(i,j))
            albmin=MIN(albmin,albedo(i,j))
        END IF
        IF(cloud_frac_vis_a(i,j) >= 0.) THEN
          vcfmax=MAX(vcfmax,cloud_frac_vis_a(i,j))
          vcfmin=MIN(vcfmin,cloud_frac_vis_a(i,j))
        END IF
      END DO
    END DO
    IF (myproc == 0) THEN
      WRITE(6,'(a)') 'Mosaicked Visible Satellite Data: '
      DO i=1, nvisfiles
        WRITE(6,'(a)') vis_fname(i)
      END DO
      WRITE(6,'(a,f10.4,a,f10.4)')                                      &
         ' min albedo:',albmin,'  max albedo:',albmax
      WRITE(6,'(a,f10.4,a,f10.4)')                                      &
         ' min vis cf:',vcfmin,'  max vis cf:',vcfmax
    END IF

  ELSE
    IF (myproc == 0) THEN
      IF (nvisfiles > 0) THEN
        WRITE(6,'(a,a,/a)')  ' ###Problem reading ',vis_fname,          &
                         ' Skipping vis processing'
      ELSE
        WRITE(6,'(a)') 'Skipping vis processing'
      END IF
    END IF
  END IF
!
  IF (cld_files == 1) THEN
    varid ='cldalb'
    varname='Cloud from Albedo'
    varunits='Fraction'
    CALL wrtvar2(nx,ny,1,cloud_frac_vis_a,varid,varname,varunits,       &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
  END IF
!
!-----------------------------------------------------------------------
!
!  Insert satellite IR data
!
!-----------------------------------------------------------------------
!
  ifield=1
  isource=1

  IF (myproc == 0) THEN
    WRITE(6,*) 'Inserting satellite IR data.'
    WRITE(6,*) 'Getting satellite IR data remapped on the model',       &
             ' grid'
  END IF
!
  CALL irmosaic(nx,ny,mx_sat,nirfiles,ir_fname,isatir,                  &
                    irtemp,tem2d2,                                      &
                    satnamir,latsatir,lonsatir,io_ir)
!
  IF (cld_files == 1) THEN
    varid='irtemp'
    varname='Cloud Top Temp'
    varunits='Degrees K'
    CALL wrtvar2(nx,ny,1, irtemp,varid,varname,varunits,                &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
  END IF
!
  IF(io_ir == 0) THEN
!    WRITE(100+mp_opt+myproc,'(a)')  ' Read successful, getting calibration coefficients'
    CALL rdircal(io_calib)
    IF(io_calib == 0) THEN
      IF (myproc == 0) WRITE(6,'(a)')  ' Read successful, calling insert_ir'
      CALL insert_ir (nx,ny,nz,latgr,longr,rland_frac,cvr_snow,         &
                    hterain,zs,z_lcl,t_3d,p_3d,                         &
                    t_sfc_k,t_gnd_k,clouds_3d,                          &
                    solar_alt,solar_ha,solar_dec,                       &
                    isatir,irtemp,cldtop_m_irt,cldtop_m,                &
                    dx,dy,sfc_sao_buffer,                               &
                    clouds_3d,tem1,lvldbg,istatus_ir)
    ELSE
      io_ir=-1
      IF (myproc == 0) WRITE(6,'(a,a,/a)')                              &
        ' !!!WARNING!!! Problem reading satellite calibration',         &
                     ' Skipping IR insertion'
    END IF
  ELSE
    IF (nirfiles > 0) THEN
      IF (myproc == 0) WRITE(6,'(a,/a)')                                &
        ' !!!WARNING!!! Problem reading satellite IR data',             &
                     ' Skipping IR insertion'
    ELSE
      IF (myproc == 0) WRITE(6,'(a)') 'Skipping IR insertion'
    END IF
  END IF
!
  IF (cld_files == 1) THEN
    IF (myproc == 0) WRITE(6,*)  'Writing cloud fraction after inserting IR data.'
    varid='cld_ir'
    varname='Clouds after IR Data'
    varunits='Fraction'
    CALL wrtvar2(nx,ny,nz,clouds_3d,varid,varname,varunits,             &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
  END IF
!
!-----------------------------------------------------------------------
!
!  Getting the 2-D maximum radar reflectivity field
!
!-----------------------------------------------------------------------
!
  DO i=1,nx
    DO j=1,ny
      dbz_max_2d(i,j) = r_missing
      DO k=1,nz
        dbz_max_2d(i,j) = MAX(dbz_max_2d(i,j),ref_mos_3d(i,j,k))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Insert radar data
!
!-----------------------------------------------------------------------
!
  330   CONTINUE
  IF (istat_radar == 1) THEN

    IF (myproc == 0) WRITE(6,*) 'Inserting radar reflectivity data.'
    CALL insert_radar(nx,ny,nz,clddiag,hterain,zs,t_3d,z_lcl,           &
                      refthr1,refthr2,hgtrefthr,                        &
                      ref_mos_3d,clouds_3d,                             &
                      cloud_base,cloud_base_buf,l_unresolved,           &
                      istatus_rad)

    IF (istatus_rad /= 1) THEN
      WRITE(6,'(1x,a)') 'Error inserting radar data. Aborting cld analysis'
      GO TO 999
    END IF
  ELSE
    IF (myproc == 0)   WRITE(6,'(1x,a)')                                &
      '!!!WARNING!!! Problem using radar data. Skipping RADAR insertion'
  END IF
!
  IF (cld_files == 1) THEN
    IF (myproc == 0) WRITE(6,*)  'Writing cloud fraction after inserting radar data.'
    varid='cldrad'
    varname='Clouds after Radar'
    varunits='Fraction'
    CALL wrtvar2(nx,ny,nz,clouds_3d,varid,varname,varunits,             &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
  END IF
!
  IF (cloudopt == 2) GO TO 340   ! Using radar data only
!
!-----------------------------------------------------------------------
!
!  Inserting the satellite visible data
!
!-----------------------------------------------------------------------
!
  IF (io_vis == 0) THEN
    IF (myproc == 0) WRITE(6,*) 'Inserting satellite VIS data.'
    CALL insert_vis (nx,ny,nz,zs,hterain,clouds_3d,                     &
                     albedo,cloud_frac_vis_a,                           &
                     vis_rad_thcvr,vis_rad_thdbz,                       &
                     istat_radar,ref_mos_3d,refthr2,dbz_max_2d,         &
                     r_missing,sfc_sao_buffer,lvldbg,istatus_vis)
!
    IF (istatus_vis /= 1) THEN
      WRITE(6,*)' Error inserting VIS data, Aborting cld analysis.'
      GO TO 999
    END IF
  ELSE
    IF (myproc == 0) WRITE(6,'(1x,a)')                                  &
      '!!!WARNING!!! Problem using satellite VIS data. Skipping VIS data insertion'
  END IF
!
  IF (cld_files == 1) THEN
    IF (myproc == 0) WRITE(6,*)  'Writing cloud fraction after inserting VIS data.'
    varid='cldvis'
    varname='Clouds after Vis'
    varunits='Fraction'
    CALL wrtvar2(nx,ny,nz,clouds_3d,varid,varname,varunits,             &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
  END IF
!
  340   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Calculate the column maximum and column total cloud cover.
!
!-----------------------------------------------------------------------
!
  DO i=1,nx
    DO j=1,ny
      DO k=1,nz
        cvr_1d (k) = clouds_3d(i,j,k)
      END DO
      CALL col_max_tot (nz,cvr_1d,cvr_max(i,j),cvr_tot(i,j)             &
                        ,istat_col)
    END DO
  END DO
  IF (cld_files == 1) THEN
    IF (myproc == 0) WRITE(6,*)  'Writing column total cloud cover.'
    varid   ='totcvr'
    varname ='Total Column Cld Cov'
    varunits='Fraction'
    CALL wrtvar2(nx,ny,1,cvr_tot,varid,varname,varunits,                &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
  END IF
!
!-----------------------------------------------------------------------
!
!  Compare the final cloud analysis with radar reflectivity
!
!-----------------------------------------------------------------------
!
  IF (istat_radar == 1 .AND. myproc == 0) THEN
    CALL compare_radar(nx,ny,nz,ref_mos_3d,dbz_max_2d                   &
                 ,cvr_max,refthr2,cloud_frac_vis_a                      &
                 ,vis_rad_thcvr,vis_rad_thdbz,r_missing)
  END IF
!
!-----------------------------------------------------------------------
!
!  Get Cloud Bases and Tops
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,*)' Calculating Cloud Ceiling, Base and Top'
!
  DO j = 1,ny-1
    DO i = 1,nx-1
      cloud_base(i,j) = default_base
      cloud_ceiling(i,j) = default_ceiling
      cloud_top(i,j)  = default_top
!
!-----------------------------------------------------------------------
!
!  Test for Cloud Base (MSL)
!
!-----------------------------------------------------------------------
!
      IF (cvr_max(i,j) >= thresh_cvr_base) THEN
        DO k = nz-1,1,-1
          IF(clouds_3d(i,j,k) < thresh_cvr_base .AND.                   &
                clouds_3d(i,j,k+1) >= thresh_cvr_base) THEN
            cloud_base(i,j) = 0.5 * (zs(i,j,k) + zs(i,j,k+1))
          END IF
        END DO ! k
      END IF ! Clouds exist in this column
!
!-----------------------------------------------------------------------
!
!  Test for Cloud Top (MSL)
!
!-----------------------------------------------------------------------
!
      IF (cvr_max(i,j) >= thresh_cvr_top) THEN
        DO k = 1,nz-1
          IF(clouds_3d(i,j,k) > thresh_cvr_top .AND.                    &
                clouds_3d(i,j,k+1) <= thresh_cvr_top) THEN
            cloud_top(i,j) = 0.5 * (zs(i,j,k) + zs(i,j,k+1))
          END IF
        END DO ! k
      END IF ! Clouds exist in this column
!
!-----------------------------------------------------------------------
!
!  Test for Cloud Ceiling (AGL)
!
!-----------------------------------------------------------------------
!
      IF (cvr_max(i,j) >= thresh_cvr_ceiling) THEN
        DO k = nz-1,1,-1
          IF(clouds_3d(i,j,k) < thresh_cvr_ceiling .AND.                &
                clouds_3d(i,j,k+1) >= thresh_cvr_ceiling) THEN
            cloud_ceiling(i,j) = 0.5 * (zs(i,j,k)+zs(i,j,k+1))          &
                                     - hterain(i,j)
          END IF
        END DO ! k
      END IF ! Clouds exist in this column

    END DO ! i
  END DO ! j
!ok 3
!
!-----------------------------------------------------------------------
!
!  Set the boundary values
!
!-----------------------------------------------------------------------
!
  IF (mp_opt == 0) THEN
    DO i = 1,nx-1
      cloud_base(i,ny) = cloud_base(i,ny-1)
      cloud_ceiling(i,ny) = cloud_ceiling(i,ny-1)
      cloud_top(i,ny)  = cloud_top(i,ny-1)
    END DO ! i
!
    DO j = 1,ny
      cloud_base(nx,j) = cloud_base(nx-1,j)
      cloud_ceiling(nx,j) = cloud_ceiling(nx-1,j)
      cloud_top(nx,j)  = cloud_top(nx-1,j)
    END DO ! i
  END IF
!
  IF (cld_files == 1) THEN
    varid='cldbas'
    varname='Cloudbase'
    varunits='Meters'
    CALL wrtvar2(nx,ny,1,cloud_base,varid,varname,varunits,             &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
    varid='cldtop'
    varname='Cloud Top Height'
    varunits='Meters'
    CALL wrtvar2(nx,ny,1,cloud_top,varid,varname,varunits,              &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
    varname='ceilng'
    varname='Cloud Ceiling'
    varunits='Meters'
    CALL wrtvar2(nx,ny,1,cloud_ceiling,varid,varname,varunits,          &
                 time,runname,dirname,hdmpfmt,hdfcompr,joindmp(FINDX_A),&
                 istatwrt)
  END IF
!
!-----------------------------------------------------------------------

  istatus=1

  999   CONTINUE
!
!-----------------------------------------------------------------------
!
  RETURN
END SUBROUTINE cloud_cv




!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE COL_MAX_TOT                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE col_max_tot (nz,cvr,cvr_max,cvr_tot,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  To find the number of cloud decks in one column of grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  06/01/96.
!
!  MODIFICATION HISTORY:
!  04/11/97  J. Zhang
!            Included adascld24.inc for ncloud
!  08/06/97  J. Zhang
!            Change adascld24.inc to adascld26.inc
!  05/05/98  J. Zhang
!            Abandon cloud grid. Using only ARPS grid.
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
!-----------------------------------------------------------------------
!
!  INPUT:
  INTEGER :: nz
!
!  INPUT:
  REAL :: cvr(nz)     ! Cloud cover in one column of grid.
!
!  OUTPUT:
  INTEGER :: istatus
  REAL :: cvr_max         ! column maximum cldcvr
  REAL :: cvr_tot         ! column integrated cldcvr
!
!  LOCAL:
  INTEGER :: nlyr         ! number of cloud decks (layers) in the column.
  INTEGER :: ilyr(nz) ! Layer index for each cloud lvl (needs)

  REAL :: a(nz)       ! max. cloud cover of each cloud layer
  REAL :: f(nz)       ! Apparnt "x-sectn" of cldlyrs seen from above
  INTEGER :: ik(nz)   ! Height level representative of cloud layers
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k,n
  REAL :: cvr_thresh,sumf
  PARAMETER (cvr_thresh = 0.1)
!
!-----------------------------------------------------------------------
!
!  Convert from cld cvr to discreet cloud layer indices (cvr to a)
!
!  Find the maximum cloud cover "a" in each cloud deck, and
!  the level index "ik" associated with the maximum cloud cover.
!
!-----------------------------------------------------------------------
!
  istatus=0
  nlyr = 0
  IF(cvr(nz) > cvr_thresh) THEN    ! cld_top at upper bndry
    nlyr = nlyr + 1
    a(nlyr) = cvr(nz)
    ik(nlyr) = nz
  END IF

  DO k = nz-1,1,-1
    IF(cvr(k) >= cvr_thresh .AND. cvr(k+1) < cvr_thresh) THEN
      nlyr = nlyr + 1        ! at top of a new cld layer
      a(nlyr) = cvr(k)
      ik(nlyr) = k
    ELSE IF(nlyr >= 1) THEN
      IF(cvr(k) > a(nlyr)) THEN
        a(nlyr) = cvr(k)   ! Max cldcvr within a layer
        ik(nlyr) = k       ! the lvl index for the max. cldcv
      END IF
    END IF

    IF(cvr(k) >= cvr_thresh) THEN      ! Still within layer
      ilyr(k) = nlyr
    ELSE                                 ! Below layer
      ilyr(k) = 0
    END IF

  END DO ! k
!
!-----------------------------------------------------------------------
!
!  Convert cloud layer fractions to "cross-section" seen from
!  satellite.  This solves for the "f" array given the "a" array
!
!-----------------------------------------------------------------------
!
  IF(nlyr == 0) THEN
    istatus=1
    cvr_max=0.0
    cvr_tot=0.0
    RETURN
  END IF

  a(1) = MIN(a(1),1.0)
  f(1) = a(1)
  sumf = f(1)
  cvr_max=a(1)
  cvr_tot=a(1)

  IF(nlyr >= 2) THEN
    DO n = 2,nlyr
      a(n) = MIN(a(n),1.0)       ! max cldcvr in one cld layer
      cvr_max = MAX(a(n),cvr_max)
      f(n) = a(n) * (1.0 - sumf) ! fraction of radiation reaches
                                 ! top of atm from each cld layer
      sumf = sumf + f(n)         ! fraction of radiation blked
                                 ! /attened by all cld lyrs above
    END DO ! n
!
!-----------------------------------------------------------------------
!
!  Calculate column integrated cloud cover as if the cloud
!  particles are uniformly (fusely) distributed in the atmosphere.
!
!-----------------------------------------------------------------------
!
    cvr_tot =sumf
!
  END IF ! nlyr
!
!-----------------------------------------------------------------------
!
  istatus=1
  RETURN
END SUBROUTINE col_max_tot
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE COMPARE_RADAR               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE compare_radar (nx,ny,nz,ref_mos_3d,dbz_max_2d,               &
                          cvr_max,ref_base,cf_vis_a,                    &
                          vis_rad_thcvr,vis_rad_thdbz,r_missing)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This routine compares the cloud and radar fields and flags
!  remaining differences that weren't caught in earlier processing
!  This routine does not alter any arrays or pass anything back,
!  it is diagnostic only.
!
!  Since it is diagnostic only, it is not necessary to be called by
!  all processes. - WYH.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:   Jian Zhang
!  07/15/96  Based on the LAPS cloud analysis code of 9/1995
!
!  MODIFICATION HISTORY:
!
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
!-----------------------------------------------------------------------
!
!  INPUT:
  INTEGER :: nx,ny,nz            ! Model grid size(dbz)
  REAL :: ref_mos_3d(nx,ny,nz) ! radar reflectivity remapped on grid
  REAL :: dbz_max_2d(nx,ny)      ! Column max. radar reflectivity (dbz)
  REAL :: cvr_max(nx,ny)         ! column maximum cloud cover
  REAL :: cf_vis_a(nx,ny)  ! sate. visible albedo derived cldcvr
  REAL :: ref_base       ! thrshld define signif. radar echo
  REAL :: vis_rad_thcvr  ! thrshld define signif. vis cldcvr (0.2)
  REAL :: vis_rad_thdbz  ! thrshld define signif. cldy radar ehco (10dbz)
  REAL :: r_missing      ! flag for missing or bad data
!
!  OUTPUT:
!  None
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
!
!-----------------------------------------------------------------------
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(a)') 'Comparing clouds and radar (_RDR)'
  WRITE(6,'(a,f5.2)') 'vis_rad_thcvr = ',vis_rad_thcvr
  WRITE(6,'(a,f5.2)') 'vis_rad_thdbz = ',vis_rad_thdbz

  DO i = 1,nx
    DO j = 1,ny
      IF(cvr_max(i,j) < vis_rad_thcvr .AND.                             &
         dbz_max_2d(i,j) > ref_base) THEN
!
!-----------------------------------------------------------------------
!
!  We have a discrepancy between the VIS and radar
!
!-----------------------------------------------------------------------
!
        IF(cvr_max(i,j) == cf_vis_a(i,j)) THEN ! CVR = VIS

          IF(dbz_max_2d(i,j) < vis_rad_thdbz) THEN
            WRITE(6,'(a,2i4,f5.2,f8.1,f5.2)')                         &
               ' VIS_RDR: cvr_mx/dbz/vis <',                          &
               i,j,cvr_max(i,j),dbz_max_2d(i,j),cf_vis_a(i,j)
          ELSE
            WRITE(6,'(a,2i4,f5.2,f8.1,f5.2)')                         &
               ' VIS_RDR: cvr_mx/dbz/vis >',                          &
               i,j,cvr_max(i,j),dbz_max_2d(i,j),cf_vis_a(i,j)
          END IF  ! dbz_max_2d < vis_rad_thdbz

        ELSE IF (cvr_max(i,j) < cf_vis_a(i,j)) THEN
!
!-----------------------------------------------------------------------
!
!  Don't blame the VIS            CVR < VIS
!
!-----------------------------------------------------------------------
!
          IF(dbz_max_2d(i,j) < vis_rad_thdbz) THEN
!
!-----------------------------------------------------------------------
!
!  We can potentially block out the radar
!
!-----------------------------------------------------------------------
!
            WRITE(6,'(a,2i4,f5.2,f8.1,f5.2)')                         &
              ' CLD_RDR: cvr_mx/dbz/vis <',                           &
              i,j,cvr_max(i,j),dbz_max_2d(i,j),cf_vis_a(i,j)
          ELSE ! Radar is too strong to block out
            WRITE(6,'(a,2i4,f5.2,f8.1,f5.2)')                         &
              ' CLD_RDR: cvr_mx/dbz/vis >',                           &
              i,j,cvr_max(i,j),dbz_max_2d(i,j),cf_vis_a(i,j)
          END IF   ! dbz_max_2d < vis_rad_thdbz

        ELSE IF (cvr_max(i,j) > cf_vis_a(i,j)) THEN
!
!-----------------------------------------------------------------------
!
!  Don't know if VIS lowered cloud cover        CVR > VIS
!  At least if difference is less than VIS "cushion"
!
!-----------------------------------------------------------------------
!
          IF(dbz_max_2d(i,j) < vis_rad_thdbz) THEN
!
!-----------------------------------------------------------------------
!
!  We can potentially block out the radar
!
!-----------------------------------------------------------------------
!
            WRITE(6,'(a,2i4,f5.2,f8.1,f8.2)')                         &
                  ' ???_RDR: cvr_mx/dbz/vis <',                       &
                  i,j,cvr_max(i,j),dbz_max_2d(i,j),cf_vis_a(i,j)
          ELSE ! Radar is too strong to block out
            WRITE(6,'(a,2i4,f5.2,f8.1,f8.2)')                         &
                  ' ???_RDR: cvr_mx/dbz/vis >',                       &
                  i,j,cvr_max(i,j),dbz_max_2d(i,j),cf_vis_a(i,j)
          END IF   ! dbz_max_2d < vis_rad_thdbz
        END IF  ! cvr_max(i,j) = cf_vis_a(i,j)
      END IF  ! cvr_max < vis_rad_thcvr & dbz_max_2d > ref_base
    END DO ! j
  END DO ! i

  RETURN
END SUBROUTINE compare_radar
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE RDIRCAL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE rdircal(io_calib)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the 10-micron IR calibration data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  03/12/2003
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER, INTENT(OUT) :: io_calib
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'adassat.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: isat,iline,iunit
  CHARACTER (LEN=80) :: chline
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0) THEN
    io_calib=-9
    CALL getunit(iunit)
    DO isat=1,nirfiles
      OPEN(iunit,file=ircalname(isat),status='old')
      DO iline=1,1001
        READ(iunit,'(a80)') chline
        IF(chline(1:1) /= '!') THEN
          READ(chline,*) centw(isat),fk1(isat),fk2(isat),               &
                                     bc1(isat),bc2(isat)
          bc2inv(isat)=1./bc2(isat)
          io_calib=0
          EXIT
        END IF
      END DO
      CLOSE(iunit)
    END DO
    CALL retunit(iunit)
  END IF

  CALL mpupdater(centw,nirfiles)
  CALL mpupdater(fk1,nirfiles)
  CALL mpupdater(fk2,nirfiles)
  CALL mpupdater(bc1,nirfiles)
  CALL mpupdater(bc2,nirfiles)
  CALL mpupdater(bc2inv,nirfiles)
  CALL mpupdatei(io_calib,1)

  RETURN
END SUBROUTINE rdircal
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE VISMOSAIC                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE vismosaic(nx,ny,mx_sat,nvisfiles,vis_fname,viscalname,       &
                     isatvis,albedo,tem2d,                              &
                     satnamvis,latsatvis,lonsatvis,io_vis)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read, calibrate and mosaic the visible satellite data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  03/12/2003
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN) :: nx,ny
  INTEGER, INTENT(IN) :: mx_sat
  INTEGER, INTENT(IN) :: nvisfiles
  CHARACTER(LEN=*), INTENT(IN) :: vis_fname(mx_sat)
  CHARACTER(LEN=*), INTENT(IN) :: viscalname(mx_sat)
  INTEGER, INTENT(OUT) :: isatvis(nx,ny)
  REAL,    INTENT(OUT) :: albedo(nx,ny)
  REAL,    INTENT(OUT) :: tem2d(nx,ny)
  CHARACTER (LEN=6), INTENT(OUT) :: satnamvis(mx_sat)
  REAL,    INTENT(OUT) :: latsatvis(mx_sat)
  REAL,    INTENT(OUT) :: lonsatvis(mx_sat)
  INTEGER, INTENT(OUT) :: io_vis
!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,isat,iline,iunit,iotem,io_cal
  INTEGER :: yrl,monl,dayl,abstsec0,abstsec,itlaunch,iday,lday,ndays
  INTEGER :: itime,isource
  REAL :: day0cal,dimrate
  REAL :: aconst,albtem
  CHARACTER (LEN=80) :: chline
  CHARACTER (LEN=10) :: ctlaunch
  CHARACTER (LEN=6)  :: varname
  INTEGER, PARAMETER :: nsecday = 86400
  INTEGER :: nxlg, nylg
  REAL, ALLOCATABLE :: tem(:,:)
  INTEGER :: istat
!
!-----------------------------------------------------------------------
!
! Include Files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  io_vis=-9
  DO j=1,ny
    DO i=1,nx
      isatvis(i,j)=1
      albedo(i,j)=-999.
    END DO
  END DO
  IF( nvisfiles == 0 ) RETURN
!
!-----------------------------------------------------------------------
!
! Initialization, current day
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
      nxlg = (nx-3) * (nproc_x) + 3
      nylg = (ny-3) * (nproc_y) + 3
      ALLOCATE(tem(nxlg,nylg),stat=istat)
  END IF

  IF (myproc == 0) THEN
    CALL ctim2abss(year,month,day,hour,minute,second, abstsec0)
    abstsec=abstsec0+nint(curtim)
    iday=abstsec/nsecday
!
!-----------------------------------------------------------------------
!
! Read first visible satellite file
!
!-----------------------------------------------------------------------
!
    IF (mp_opt > 0 ) THEN
      CALL rdsatfld(nxlg,nylg,1,vis_fname(1),                           &
                  satnamvis(1),latsatvis(1),lonsatvis(1),               &
                  itime,isource,varname,tem,io_vis)
    ELSE
      CALL rdsatfld(nx,ny,1,vis_fname(1),                               &
                  satnamvis(1),latsatvis(1),lonsatvis(1),               &
                  itime,isource,varname,albedo,io_vis)
    END IF
  END IF
  IF (mp_opt > 0) CALL mpisplit2d(tem,nx,ny,albedo)
  CALL mpupdatei(io_vis,1)

  !
  !-----------------------------------------------------------------------
  !
  ! Read calibration file
  !
  !-----------------------------------------------------------------------
  !
  IF (myproc == 0) THEN
    CALL getunit(iunit)
    OPEN(iunit,file=viscalname(1),status='old')
    DO iline=1,1001
      READ(iunit,'(a80)') chline
      IF(chline(1:1) /= '!') THEN
        READ(chline,*) ctlaunch,day0cal,dimrate
        EXIT
      END IF
    END DO
    CLOSE(iunit)
    CALL retunit(iunit)
  !
  !-----------------------------------------------------------------------
  !
  ! Calculate calibration constant
  !
  !-----------------------------------------------------------------------
  !
    READ(ctlaunch,'(i4,1x,i2,1x,i2)') yrl,monl,dayl
    CALL ctim2abss(yrl,monl,dayl,12,0,0,itlaunch)
    lday=itlaunch/nsecday
    ndays=max(0,(iday-lday))
    aconst=day0cal*(1.0+(ndays*dimrate))
  END IF

  CALL mpupdater(aconst,1)

  !
  !-----------------------------------------------------------------------
  !
  ! Apply calibration
  !
  !-----------------------------------------------------------------------
  !
  DO j=1,ny
    DO i=1,nx
      IF ( albedo(i,j) > 0. ) THEN
        albedo(i,j)=min((albedo(i,j)*aconst),1.0)
      END IF
    END DO
  END DO

  IF( nvisfiles == 1 ) THEN
    IF (myproc == 0 .AND. mp_opt > 0) DEALLOCATE(tem)
    RETURN
  END IF
!
  DO isat=2,nvisfiles
!
!-----------------------------------------------------------------------
!
! Read next visible satellite file
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0) THEN
      IF (mp_opt > 0) THEN
        CALL rdsatfld(nxlg,nylg,1,vis_fname(isat),satnamvis(isat),      &
                latsatvis(isat),lonsatvis(isat),                        &
                itime,isource,varname,tem,iotem)
      ELSE
        CALL rdsatfld(nx,ny,1,vis_fname(isat),satnamvis(isat),          &
                latsatvis(isat),lonsatvis(isat),                        &
                itime,isource,varname,tem2d,iotem)
      END IF
    END IF
    IF (mp_opt > 0) CALL mpisplit2d(tem,nx,ny,tem2d)
    CALL mpupdatei(iotem,1)
    IF( iotem == 0 ) THEN
!
!-----------------------------------------------------------------------
!
! Read calibration file
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0) THEN
        CALL getunit(iunit)
        OPEN(iunit,file=viscalname(isat),status='old')
        DO iline=1,1001
          READ(iunit,'(a80)') chline
          IF(chline(1:1) /= '!') THEN
            READ(chline,*) ctlaunch,day0cal,dimrate
            EXIT
          END IF
        END DO
        CLOSE(iunit)
        CALL retunit(iunit)
!
!-----------------------------------------------------------------------
!
! Calculate calibration constant
!
!-----------------------------------------------------------------------
!
        READ(ctlaunch,'(i4,1x,i2,1x,i2)') yrl,monl,dayl
        CALL ctim2abss(yrl,monl,dayl,12,0,0,itlaunch)
        lday=itlaunch/nsecday
        ndays=max(0,(iday-lday))
        aconst=day0cal*(1.0+(ndays*dimrate))
      END IF
      CALL mpupdater(aconst,1)
!
!-----------------------------------------------------------------------
!
! Apply calibration
!
!-----------------------------------------------------------------------
!
      DO j=1,ny
        DO i=1,nx
          IF(tem2d(i,j) >= 0. ) THEN
             albtem=max((tem2d(i,j)*aconst),1.0)
             IF(albtem > albedo(i,j)) THEN
                albedo(i,j)=albtem
                isatvis(i,j)=isat
             END IF
          END IF
        END DO
      END DO
    END IF
  END DO
  IF (myproc == 0 .AND. mp_opt > 0) DEALLOCATE(tem)

  RETURN
END SUBROUTINE vismosaic
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE IRMOSAIC                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE irmosaic(nx,ny,mx_sat,nirfiles,ir_fname,isatir,              &
                    irtemp,tem2d2,                                      &
                    satnamir,latsatir,lonsatir,io_ir)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read and mosaic the IR (10.7 micron) satellite data files.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  03/12/2003
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN) :: nx,ny
  INTEGER, INTENT(IN) :: mx_sat
  INTEGER, INTENT(IN) :: nirfiles
  CHARACTER (LEN=256), INTENT(IN) :: ir_fname(mx_sat)
  INTEGER, INTENT(OUT) :: isatir(nx,ny)
  REAL,    INTENT(OUT) :: irtemp(nx,ny,2)
  REAL,    INTENT(OUT) :: tem2d2(nx,ny,2)
  CHARACTER (LEN=6), INTENT(OUT) :: satnamir(mx_sat)
  REAL,    INTENT(OUT) :: latsatir(mx_sat)
  REAL,    INTENT(OUT) :: lonsatir(mx_sat)
  INTEGER, INTENT(OUT) :: io_ir
!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,isat,iotem
  INTEGER :: itime,isource
  CHARACTER (LEN=6) :: varname
  INTEGER :: nxlg, nylg
  REAL, ALLOCATABLE :: tem(:,:,:)
  INTEGER :: istat
!
!-----------------------------------------------------------------------
!
! Initialize output data arrays
!
!-----------------------------------------------------------------------
!
  io_ir=-9
  isatir=1
  irtemp=-999.

  IF(nirfiles == 0) RETURN
!
!-----------------------------------------------------------------------
!
! Read the first IR satellite file
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    nxlg = (nx-3) * (nproc_x) + 3
    nylg = (ny-3) * (nproc_y) + 3
    ALLOCATE(tem(nxlg,nylg,2),stat=istat)
  END IF

  IF (myproc == 0) THEN
    IF (mp_opt > 0) THEN
      CALL rdsatfld(nxlg,nylg,2,ir_fname(1),satnamir(1),                &
                latsatir(1),lonsatir(1),                                &
                itime,isource,varname,                                  &
                tem,io_ir)
    ELSE
      CALL rdsatfld(nx,ny,2,ir_fname(1),satnamir(1),                    &
                latsatir(1),lonsatir(1),                                &
                itime,isource,varname,                                  &
                irtemp,io_ir)
    ENDIF
  END IF
  CALL mpupdatei(io_ir,1)
  IF (mp_opt > 0) THEN
!    CALL mpisplit2d(tem(:,:,1),nx,ny,irtemp(1,1,1))
!    CALL mpisplit2d(tem(:,:,2),nx,ny,irtemp(1,1,2))
    CALL mpisplit3d(tem,nx,ny,2,irtemp)
  END IF

  IF(nirfiles == 1) THEN
    IF (myproc == 0 .AND. mp_opt > 0) DEALLOCATE(tem)
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
! Read and merge second and subsequent files.
!
!-----------------------------------------------------------------------
!
  DO isat=2,nirfiles
    IF (myproc == 0) THEN
      IF (mp_opt > 0) THEN
        CALL rdsatfld(nxlg,nylg,2,ir_fname(isat),satnamir(isat),        &
                latsatir(isat),lonsatir(isat),                          &
                itime,isource,varname,                                  &
                tem,iotem)
      ELSE
        CALL rdsatfld(nx,ny,2,ir_fname(isat),satnamir(isat),            &
                latsatir(isat),lonsatir(isat),                          &
                itime,isource,varname,                                  &
                tem2d2,iotem)
      END IF
    END IF
    IF (mp_opt > 0) THEN
!      CALL mpisplit2d(tem(:,:,1),nx,ny,tem2d2(1,1,1))
!      CALL mpisplit2d(tem(:,:,2),nx,ny,tem2d2(1,1,2))
       CALL mpisplit3d(tem,nx,ny,2,tem2d2)
    END IF
    IF( io_ir == 0 ) THEN
      DO j=1,ny
        DO i=1,nx
          IF(irtemp(i,j,1) < 0.0) THEN
            isatir(i,j)=isat
            irtemp(i,j,1)=tem2d2(i,j,1)
            irtemp(i,j,2)=tem2d2(i,j,2)
          ELSE IF( tem2d2(i,j,1) > 0.0 .and.                            &
                  (tem2d2(i,j,1) < irtemp(i,j,1)) ) THEN
            isatir(i,j)=isat
            irtemp(i,j,1)=tem2d2(i,j,1)
            irtemp(i,j,2)=tem2d2(i,j,2)
          END IF
        END DO
      END DO
    END IF
  END DO
  IF (myproc == 0 .AND. mp_opt > 0) DEALLOCATE(tem)
  RETURN
END SUBROUTINE irmosaic
