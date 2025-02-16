!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE INSERT_SAO1                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE insert_sao1(nx,ny,nz,dx,dy,xs,ys,zs,hterain,t_sfc_k,         &
           nxlg,nylg,xslg,yslg,                                         &
           cldcv,wtcldcv,t_modelfg,cf_modelfg,                          &
           nobsng,indexsng,stn,isrcsng,obstype,xsng,ysng,ista_snd,      &
           obstime,latsta,lonsta,elevsta,                               &
           kcloud,store_amt,store_hgt,                                  &
           l_stn_clouds,n_cld_snd,cld_snd,wt_snd,                       &
           i_snd,j_snd,                                                 &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Ingest the SAO cloud coverage observations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  02/96.  Based on the LAPS cloud analysis code by Steve Albers,
!          1995.
!
!  MODIFICATION HISTORY:
!
!  02/06/96  J. Zhang
!            Modified for ADAS grid. Added documents.
!  03/14/97  J. Zhang
!            Clean up the code and implemented for the official
!            arps4.2.4 version
!  04/15/97  J. Zhang
!            Added error message output for the case when the
!            actual # of cloud soundings (N_CLD_SND) exceeds the
!            maximum allowed # (MAX_CLD_SND) defined in
!            adascld24.inc.
!  08/06/97  J. Zhang
!            Change adascld24.inc to adascld25.inc
!  09/09/97  J. Zhang
!            Assign a weight to a SAO cloud sounding base on
!            the cloud coverage amount.
!  09/11/97  J. Zhang
!            Change adascld25.inc to adascld26.inc
!  02/17/98  Keith Brewster
!            Changed logic to skip cloud layers with below zero
!            height -- missing, in other words.
!  05/05/98  J. Zhang
!            Abandoned the cloud grid, using the ARPS grid instead.
!  05/01/98  Keith Brewster
!            Removed abort for missing data, instead data are
!            checked for missing at every step.
!  04/10/03  K. Brewster
!            Modified some variable names for ADAS consistency.
!
!-----------------------------------------------------------------------
!
!  INCLUDE: (from adas.inc)
!
!  mx_sng            ! max. possible # of surface obs.
!  max_cld_snd       ! max. possible # of stations
                     ! with cloud coverage reports.
!
!  INPUT:
!
!  nx,ny,nz          ! Number of ADAS grid pts.in 3 directions.
!
!  l_stn_clouds       ! Using SAO stations' cloud obs?
!
!c ARPS grid variables.
!
!   xs     (nx)          ! The x-coord. of the physical and
                         ! computational grid. Defined at p-point.
!   ys     (ny)          ! The y-coord. of the physical and
                         ! computational grid. Defined at p-point.
!   zs    (nx,ny,nz)     ! The physical height coordinate defined at
                         ! p-point of the staggered grid.
!   hterain (nx,ny)      ! The height of the terrain (equivalent
!
!   First guess fields
!
!   cf_modelfg (nx,ny,nz)  ! first guess cloud cover field
!   t_modelfg (nx,ny,nz)   ! first guess temperature field
!   t_sfc_k (nx,ny)            ! surface temperature field
!
!  cldcv (nx,ny,nz)     3D gridded fractional cloud cover analysis.
!                           (when input, it's the first guess field)
!  wtcldcv (nx,ny,nz)   weights assigned to gridded fractional
!                           cloud cover analysis. (when input, it's
!                           the weights for first guess field)
!
!   Single-level (e.g., surface) station variables
!
!
!   stn (mx_sng)         ! station name of single-lvel data
!   obstype (mx_sng)     ! names of sfc sources
!   isrcsng (mx_sng)     ! whether stations are used
!   obstime (mx_sng)     ! time for the observation.
!   latsta (mx_sng)      ! latitude of single-level data
!   lonsta (mx_sng)      ! longitude of single-level data
!   elevsta (mx_sng)     ! height of single-level data
!   xsng (mx_sng)        ! x location of single-level data
!   ysng (mx_sng)        ! y location of single-level data
!
!   kcloud (mx_sng)      ! number of obs. cloud layers.
!   store_amt (mx_sng,5) ! cloud coverage (ea. layer,
                         ! ea. station).
!   store_hgt (mx_sng,5) ! height of obs. cloud layers.
!
!  OUTPUT:
!
!  istatus            The flag indicating process status.
!
!  n_cld_snd         number of cloud soundings created
!  cld_snd (max_cld_snd,nz)   Cld snding obtained from SAO data.
!  wt_snd (max_cld_snd,nz)    weights assigned to cloud sounding
!                                obtained from SAO data.
!  i_snd (max_cld_snd)       i-location of each cloud sounding
!                                  station in ADAS grid
!  j_snd (max_cld_snd)       j-location of each cloud sounding
!                                  station in ADAS grid
!
!  LOCAL:
!
!  name_array (max_cld_snd)  station names (first letter) for
!                                  each cloud sounding
!  ista_snd (max_cld_snd)    index of cloud sounding stations
!  cvr_snd (max_cld_snd)     column integrated cloud cover
!  cloud_top(nx,ny)    Analyzed cloud top heights (m ASL).
!  cloud_base(nx,ny)   Analyzed cloud base heights (m ASL).
!  cloud_ceiling(nx,ny)   Analyzed cloud ceiling heights (m ASL).
!
!-----------------------------------------------------------------------
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
  INCLUDE 'adas.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  INCLUDE: (from adas.inc)
!
!  integer mx_sng            ! max. possible # of surface obs.
!  integer max_cld_snd       ! max. possible # of stations
                             ! with cloud coverage reports.
!
!  INPUT:

  INTEGER :: nx,ny,nz
  INTEGER :: nxlg,nylg

  REAL :: dx,dy              ! ADAS grid spacing
  LOGICAL :: l_stn_clouds    ! Using SAO stns' cloud obs? (or bkgrnd)

  REAL :: xs    (nx)         ! The x-coord. of the physical and
                             ! computational grid. Defined at p-point.
  REAL :: ys    (ny)         ! The y-coord. of the physical and
                             ! computational grid. Defined at p-point.
  REAL :: zs    (nx,ny,nz)   ! The physical height coordinate defined
                             ! at p-point of the staggered grid.
  REAL :: xslg(nxlg)
  REAL :: yslg(nylg)
  REAL :: hterain (nx,ny)    ! The height of the terrain (equivalent
!
!  First guess fields
!
  REAL :: cf_modelfg (nx,ny,nz)  ! first guess cloud cover field
  REAL :: t_modelfg (nx,ny,nz)   ! first guess temperature field
  REAL :: t_sfc_k (nx,ny)            ! surface temperature field
!
  REAL :: cldcv (nx,ny,nz)      ! 3D gridded frac cld cv analysis.
                                  ! (it's also bkgrnd field when input)
  REAL :: wtcldcv (nx,ny,nz)    ! wts assgned to cld cvr analysis
                                  ! (it's also bkgrnd field when input)
!
!  Surface cloud coverage reports
!
  INTEGER :: nobsng
  INTEGER :: indexsng(nobsng)    ! "owner" of the obs
  INTEGER :: indexsng_wrk(nobsng)! work array
  CHARACTER (LEN=5) :: stn (mx_sng)    ! station name of single-lvel data
  INTEGER :: isrcsng (mx_sng)    ! is station used?
  CHARACTER (LEN=8) :: obstype (mx_sng)! names of sfc sources
  INTEGER :: obstime (mx_sng)    ! time for the observation.
  REAL :: xsng (mx_sng)          ! x location of single-level data
  REAL :: ysng (mx_sng)          ! y location of single-level data
  REAL :: latsta (mx_sng)        ! latitude of single-level data
  REAL :: lonsta (mx_sng)        ! longitude of single-level data
  REAL :: elevsta (mx_sng)       ! height of single-level data
!
  INTEGER :: kcloud (mx_sng)     ! number of cloud layers.
  CHARACTER (LEN=4) :: store_amt (mx_sng,5) ! cld coverage (ea.lyr, ea. stn).
  REAL :: store_hgt (mx_sng,5)        ! height of cloud layers.
!
!  OUTPUT:
!
  INTEGER :: istatus                ! Flag for process status
!
  INTEGER :: n_cld_snd              ! # of cloud snds created
  REAL :: cld_snd (max_cld_snd,nz)   ! Cld snd obtained from SAO data.
  REAL :: wt_snd (max_cld_snd,nz)    ! wgt for SAO cld snd
  INTEGER :: i_snd (max_cld_snd)    ! i-lctn of cld snd stn in ADAS grid
  INTEGER :: j_snd (max_cld_snd)    ! j-lctn of cld snd stn in ADAS grid
!
!  LOCAL:
!
  CHARACTER (LEN=1) :: name_array (nx,ny) ! Cld snd stn name (1st letter)
  INTEGER :: ista_snd (max_cld_snd) ! index of cld snd stns
  REAL :: cvr_snd (max_cld_snd)     ! column integrated cloud cover
!
!-----------------------------------------------------------------------
!
!  Empirical factors (parameters) for cloud cover analysis
!
!-----------------------------------------------------------------------
!
!c Using SAO stations slightly outside the ARPS domain.

  INTEGER :: ix_low,ix_high,iy_low,iy_high
  INTEGER ::  ix_low_lg,ix_high_lg,iy_low_lg,iy_high_lg
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL :: ht_defined           ! the reachable height of the SAOs
  REAL :: ht_base,ht_top       ! base and top heights of a cloud layer
  INTEGER :: n_analyzed        ! # of stations analyzed in this routine
  INTEGER :: i,k,l
  REAL :: ri,rj                ! i-,j-lctns of each stn in ADAS grid
  REAL :: rilg, rjlg           ! large grid version
  INTEGER :: iloc,jloc         ! i- and j-index of each sta.
  INTEGER :: iloclg,jloclg     ! large grid version
  INTEGER :: igrd,jgrd         ! i- and j-index of each sta. in grid
  LOGICAL :: keep              !  do we keep the data point?
  REAL :: cover                ! cloud fractional cover values
  REAL :: cld_thk

  INTEGER :: k_ceil

  LOGICAL :: l_out_of_bounds,l_dry,cldprt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0
  IF (nobsng < 1) THEN
    IF (myproc == 0)                                                        &
      WRITE(6,'(a)') '## No SAO data available. Returning from INSERT_SAO...'
    istatus = 1
    RETURN
  END IF
  IF (myproc == 0)                                                          &
    WRITE(6,'(a,i8,a)') ' Inserting SAO data from ',nobsng,' stations.'
  cldprt=(clddiag == 1)
!
  ix_low  = 1  - i_perimeter
  ix_high = nx + i_perimeter
  iy_low  = 1  - i_perimeter
  iy_high = ny + i_perimeter

!
!-----------------------------------------------------------------------
!
!  If we're MPI, everyone needs to be on the same page.  We have to use
!  large domain values.
!
!-----------------------------------------------------------------------
!

  IF (mp_opt > 0 ) THEN
    ix_low_lg = ix_low
    ix_high_lg = nxlg + i_perimeter
    iy_low_lg = iy_low
    iy_high_lg = nylg + i_perimeter
  ENDIF

!
!-----------------------------------------------------------------------
!
!  Initialize the cloud sounding array.
!
!-----------------------------------------------------------------------
!
  DO i=1,max_cld_snd
    DO k=1,nz
      cld_snd(i,k) = 0.0
      wt_snd(i,k) = 0.0
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Initially assign a unreal number for the sensor's reachable
!  height.
!
!-----------------------------------------------------------------------
!
  ht_defined = 99999.
  n_analyzed = 0
!
!-----------------------------------------------------------------------
!
!  Loop through the stations
!
!-----------------------------------------------------------------------
!
  DO i=1,nobsng

!
!-----------------------------------------------------------------------
!
!   Has the station been marked not to use?  Ob may be bad or combined
!   via "suprob".  These have been flagged with "isrcsng=-1".
!
!-----------------------------------------------------------------------
!


    IF( obstype(i)(1:4) == 'meso') GO TO 125
    IF( obstype(i)(1:4) == 'armmn') GO TO 125
    IF( obstype(i)(1:4) == 'isws') GO TO 125
    IF( obstype(i)(1:4) == 'coagmet') GO TO 125
    IF( obstype(i)(1:4) == 'hplains') GO TO 125
    IF( obstype(i)(1:4) == 'wpdn') GO TO 125
    IF( isrcsng(i) == -1) GO TO 125

!
!-----------------------------------------------------------------------
!
!   Place station at proper ADAS grid point.
!
!-----------------------------------------------------------------------
!
    ri = 1. + (xsng(i) - xs(1))/dx
    rj = 1. + (ysng(i) - ys(1))/dy
    iloc = nint(ri)
    jloc = nint(rj)

!
!   If we're MPI, check and see if the ob would be used in the domain,
!   if not, reject.  This preserves the bookkeeping.
!

    keep = .TRUE.
    IF (mp_opt > 0 ) THEN
      rilg = 1. + (xsng(i) - xslg(1))/dx
      rjlg = 1. + (ysng(i) - yslg(1))/dy
      iloclg = nint(rilg)
      jloclg = nint(rjlg)
      IF( iloclg < ix_low_lg .OR. iloclg > ix_high_lg                     &
          .OR. jloclg < iy_low_lg .OR. jloclg > iy_high_lg) THEN
        keep = .FALSE.
      END IF
    ELSE
      IF( iloc < ix_low .OR. iloc > ix_high                               &
          .OR. jloc < iy_low .OR. jloc > iy_high) THEN
        keep = .FALSE.
      END IF
    END IF

    IF ( .NOT. keep ) THEN
!      IF (myproc == 0)                                                   &
!      write(6,*) 'Station is outside domain ',stn(i)
       GO TO 125
    END IF

    IF( kcloud(i) == 0) THEN        ! Kick out AMOS stations not
                                    ! reporting clouds this time
!      write(6,*)' No cloud layers - probably AMOS -',
!    :    ' goto end of loop ',stn(i),obstype(i)
      GO TO 125
    END IF

    n_analyzed = n_analyzed + 1
    n_cld_snd = n_cld_snd + 1
    IF(n_cld_snd > max_cld_snd) THEN
      PRINT*,'##ERROR: actual # of cldsnd=',n_cld_snd,                  &
             ' exceeds max. # allowed =',max_cld_snd
      PRINT*,'Please increase MAX_CLD_SND in the .inc file.'
      PRINT*,'ABORTING from INSERT_SAO......'
      STOP
    END IF

    IF (mp_opt > 0) THEN
      indexsng_wrk(n_cld_snd) = indexsng(i)
!
!    Only do the responsibility computations if we're really "on the edge",
!    as this is what non-MPI does for points outside the grid.
!

      IF (loc_x == 1 .OR. loc_x == nproc_x .OR.                         &
       loc_y == 1 .OR. loc_y == nproc_y) THEN
!
!   If we're outside the grid, someone has to take responsibility for us,
!   otherwise, MPI and non-MPI will generate different solutions!
!
!   Only one processor should take responsibility, and we need to make sure
!   it is an edge processor!
!

!  Top or bottom edge

        IF (loc_y == 1 .OR. loc_y == nproc_y) THEN
           IF (loc_x == 1) THEN
             ix_high = nx
           ELSE IF (loc_x == nproc_x) THEN
             ix_low = 1
           ELSE
             ix_low = 1
             ix_high = nx
           END IF
        END IF

!  Left or right edge

        IF (loc_x == 1 .OR. loc_x == nproc_x) THEN
           IF (loc_y == 1) THEN
             iy_high = ny
           ELSE IF (loc_y == nproc_y) THEN
             iy_low = 1
           ELSE
             iy_low = 1
             iy_high = ny
           END IF
        END IF

        IF (indexsng_wrk(n_cld_snd) < 0) THEN
          IF( iloc >= ix_low .AND. iloc <= ix_high                      &
          .AND. jloc >= iy_low .AND. jloc <= iy_high) THEN
!       write(6,*) 'Processor ',myproc,' is responsible for ob ',n_cld_snd,' ',stn(i)
!       call flush(6)
            indexsng_wrk(n_cld_snd) = myproc
          END IF
        END IF
      END IF
    END IF

!
!-----------------------------------------------------------------------
!
!   For the mpi case...even if we don't "own" a station, there are still
!   computations that have to be made.  Specifically, "i_snd" and "j_snd"
!   have to be done
!
!-----------------------------------------------------------------------
!

    ista_snd(n_cld_snd) = i

    IF(obstype(i)(5:5) /= ' ' .AND. obstype(i)(4:7) /= 'AMOS') THEN
!
!-----------------------------------------------------------------------
!
!      Automated Station (12000' limit)
!
!-----------------------------------------------------------------------
!
      ht_defined = elevsta(i) + 12000./3.281
    ELSE
      ht_defined = 99999.
    END IF

    IF(cldprt) THEN
      WRITE(6,1,ERR=110) stn(i),latsta(i),lonsta(i)                     &
                       ,iloc,jloc,kcloud(i)                             &
                       ,obstype(i),ht_defined
      1         FORMAT(1X,a5,' lat=',f7.2,' lon=',f7.2,' i=',i3,' j=',i3 &
               ,' kld=',i1,' obsty=',a8,' h_def=',f8.0)
      110       WRITE(6,2,ERR=3)                                        &
               (store_amt(i,l),store_hgt(i,l),l=1,kcloud(i))
      2         FORMAT(1X,5(a4,f8.0))
      3         CONTINUE
    END IF

    IF( iloc < 1 .OR. iloc > (nx-1) .OR. jloc < 1 .OR. jloc > (ny-1) ) THEN
      l_out_of_bounds = .true.
      igrd=MIN(MAX(iloc,1),(nx-1))
      jgrd=MIN(MAX(jloc,1),(ny-1))
    ELSE
      l_out_of_bounds = .false.
      igrd=iloc
      jgrd=jloc
      name_array(iloc,jloc            )=stn(i)(1:1)
      name_array(iloc,MIN(jloc+1,ny))  =stn(i)(1:1)
      name_array(iloc,MAX(jloc-1,1   ))=stn(i)(1:1)
    END IF

    cvr_snd (n_cld_snd) = 0.          ! column intg. cloud cover.
!
!-----------------------------------------------------------------------
!
!  For each station, loop through each cloud layer observed
!  by the SAO.
!
!-----------------------------------------------------------------------
!
    DO l = 1,kcloud(i)

      cover = 0.

      IF (store_hgt(i,l) < 0.0 .AND. store_amt(i,l) /= ' CLR' ) THEN
!        ht_base = cld_base(l)
!        print*,' stn',i,' level',l,' cld.ht.=',store_hgt(i,l)
!    :            ,'  adj. cld. ht.=',ht_base
!        ht_base = elevsta(i) + ht_base
        CYCLE
      ELSE
        ht_base = elevsta(i) + store_hgt(i,l)
      END IF

      IF( ht_base > ht_defined+1.) THEN

        IF( store_amt(i,l) /= ' CLR') THEN  ! Clouds

          WRITE(6,*)' Error, inconsistent SAO data, cloud base is'      &
              //' reported to be too high for this sensor type'
          WRITE(6,*) ht_base,ht_defined,obstype(i)

          WRITE(6,*)' Please check cloud layer heights in the LSO'      &
              //' file to see that they are compatable with the'
          WRITE(6,*)' types of sensors used.'
          istatus = 0
          RETURN

        ELSE                             ! CLR
          WRITE(6,*)' Note, CLR sky cloud base does not'                &
              //' reflect this sensor type'
          WRITE(6,*) ht_base,ht_defined,obstype(i)

        END IF                            ! Clouds

      END IF          !"ht_base > ht_defined+1"
!
!-----------------------------------------------------------------------
!
!  CLOUDS ARE NOW IN MSL
!  Fill in clear for entire column for SAO or up to ht_base for AWOS
!
!-----------------------------------------------------------------------
!
      IF( store_amt(i,l) == ' CLR') THEN
        cover=.01
        DO k=1,nz

          IF( zs(igrd,jgrd,k) <= ht_base                                &
                .AND. zs(igrd,jgrd,k) <= ht_defined ) THEN

            CALL spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd           &
                         ,max_cld_snd,nz,iloc,jloc,k,cover,1.)
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Assign cloud cover values for entire column for SAO
!  or up to ht_base for AWOS. First check if there is dry
!  layer or inversion.
!
!-----------------------------------------------------------------------
!
      IF( store_amt(i,l) == '-SCT' ) THEN
        cover=.15 ! .2
        ht_top=ht_base+cld_thk(ht_base-elevsta(i),cover)

        DO k=2,nz-2

          IF(zs(igrd,jgrd,k) >= ht_base .AND. zs(igrd,jgrd,k) <= ht_top) THEN
!
!-----------------------------------------------------------------------
!
!            Check if inversion or dry layer exist.
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,0,l_dry)

            IF(.NOT. l_dry) THEN
              CALL spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd         &
                    ,max_cld_snd,nz,iloc,jloc,k,cover,1.)
            END IF

          ELSE IF(zs(igrd,jgrd,k) < ht_base) THEN
!
!-----------------------------------------------------------------------
!
!            Initialize the modify sounding routine
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain,          &
                   t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs,                     &
                   ht_base,ht_top,1,l_dry)

          END IF  ! ht_base < zs < ht_top
        END DO    ! k = 1, nz
      END IF

      IF( store_amt(i,l) == ' SCT') THEN
        cover=.25 ! .3
        ht_top=ht_base+cld_thk(ht_base-elevsta(i),cover)

        DO k=2,nz-2

          IF(zs(igrd,jgrd,k) >= ht_base .AND. zs(igrd,jgrd,k) <= ht_top) THEN
!
!-----------------------------------------------------------------------
!
!            Check if inversion or dry layer exist.
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,0,l_dry)

            IF(.NOT. l_dry) THEN
              CALL spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd         &
                    ,max_cld_snd,nz,iloc,jloc,k,cover,1.)
            END IF

          ELSE IF(zs(igrd,jgrd,k) < ht_base) THEN
!
!-----------------------------------------------------------------------
!
!            Initialize the modify sounding routine
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,1,l_dry)

          END IF
        END DO
      END IF

      IF( store_amt(i,l) == '-BKN' ) THEN
        cover=.4 ! .5
        ht_top=ht_base+cld_thk(ht_base-elevsta(i),cover)
        DO k=2,nz-2

          IF(zs(igrd,jgrd,k) >= ht_base .AND. zs(igrd,jgrd,k) <= ht_top) THEN
!
!-----------------------------------------------------------------------
!
!            Check if inversion or dry layer exist.
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,0,l_dry)

            IF(.NOT. l_dry)THEN
              CALL spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd         &
                    ,max_cld_snd,nz,iloc,jloc,k,cover,1.)
            END IF

          ELSE IF(zs(igrd,jgrd,k) < ht_base) THEN
!
!-----------------------------------------------------------------------
!
!            Initialize the modify sounding routine
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,1,l_dry)

          END IF
        END DO
      END IF

      IF( store_amt(i,l) == ' BKN' ) THEN
        cover=.7
        ht_top = ht_base + cld_thk(ht_base-elevsta(i),cover) ! 1500.

        DO k=2,nz-2

          IF(zs(igrd,jgrd,k) >= ht_base .AND. zs(igrd,jgrd,k) <= ht_top) THEN
!
!-----------------------------------------------------------------------
!
!            Check if inversion or dry layer exist.
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,0,l_dry)

            IF(.NOT. l_dry) THEN
              CALL spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd         &
                    ,max_cld_snd,nz,iloc,jloc,k,cover,1.)
            END IF

          ELSE IF(zs(igrd,jgrd,k) < ht_base) THEN
!
!-----------------------------------------------------------------------
!
!            Initialize the modify sounding routine
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,1,l_dry)

          END IF
        END DO
      END IF

      IF( store_amt(i,l) == '-OVC') THEN
        cover=.6 ! .9
        ht_top=ht_base+cld_thk(ht_base-elevsta(i),cover)
        DO k=1,nz

          IF(zs(igrd,jgrd,k) >= ht_base .AND. zs(igrd,jgrd,k) <= ht_top) THEN
!
!-----------------------------------------------------------------------
!
!            Check if inversion or dry layer exist.
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,0,l_dry)

            IF(.NOT. l_dry) THEN
              CALL spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd         &
                    ,max_cld_snd,nz,iloc,jloc,k,cover,1.)
            END IF

          ELSE IF(zs(igrd,jgrd,k) < ht_base) THEN
!
!-----------------------------------------------------------------------
!
!            Initialize the modify sounding routine
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,1,l_dry)

          END IF
        END DO
      END IF

      IF( store_amt(i,l) == ' OVC') THEN
        cover=1.01
        ht_top = ht_base + cld_thk(ht_base-elevsta(i),cover) ! 1500.

        DO k=2,nz-2

          IF(zs(igrd,jgrd,k) >= ht_base .AND. zs(igrd,jgrd,k) <= ht_top) THEN
!
!-----------------------------------------------------------------------
!
!            Check if inversion or dry layer exist.
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,0,l_dry)

            IF(.NOT. l_dry)THEN
              CALL spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd         &
                    ,max_cld_snd,nz,iloc,jloc,k,cover,1.)
            END IF

          ELSE IF(zs(igrd,jgrd,k) < ht_base) THEN
!
!-----------------------------------------------------------------------
!
!            Initialize the modify sounding routine
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,1,l_dry)

          END IF
        END DO
      END IF

      IF( store_amt(i,l) == ' X  ') THEN
        cover=1.01
        ht_top = ht_base + cld_thk(ht_base-elevsta(i),cover) ! 1500.

        DO k=2,nz-2

          IF(zs(igrd,jgrd,k) >= ht_base .AND. zs(igrd,jgrd,k) <= ht_top) THEN
!
!-----------------------------------------------------------------------
!
!            Check if inversion or dry layer exist.
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,0,l_dry)

            IF(.NOT. l_dry) THEN
              CALL spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd         &
                    ,max_cld_snd,nz,iloc,jloc,k,cover,1.)
            END IF
          ELSE IF(zs(igrd,jgrd,k) < ht_base) THEN
!
!-----------------------------------------------------------------------
!
!            Initialize the modify sounding routine
!
!-----------------------------------------------------------------------
!
            CALL modify_sounding(cf_modelfg,t_modelfg,hterain           &
                  ,t_sfc_k,igrd,jgrd,k,nx,ny,nz,zs                      &
                  ,ht_base,ht_top,1,l_dry)

          END IF
        END DO
      END IF

      cvr_snd(n_cld_snd) = 1. - ((1. - cvr_snd(n_cld_snd)) * cover)

    END DO
!
!-----------------------------------------------------------------------
!
!    Locate the highest ceiling
!
!-----------------------------------------------------------------------

    k_ceil = nz-2
    IF( l_stn_clouds) THEN
      DO k=nz-2,2,-1
        IF( wt_snd(n_cld_snd,k) == 1.00 .AND. cld_snd(n_cld_snd,k) > 0.5) THEN
          k_ceil = k
          GO TO 1001
        END IF
      END DO
    ELSE
      DO k=nz-2,2,-1
        IF( wtcldcv(igrd,jgrd,k) == 1.00 .AND. cldcv(igrd,jgrd,k) > 0.5 ) THEN
          k_ceil = k
          GO TO 1001
        END IF
      END DO
    END IF
!
!-----------------------------------------------------------------------
!
!  Fill in other clear layers outside of clouds, below the ceiling,
!  and within defined height range of sensor.
!
!-----------------------------------------------------------------------
!
    1001    cover = .01

    IF( l_stn_clouds) THEN
      DO k=2,k_ceil
        IF( wt_snd(n_cld_snd,k) /= 1.00                                 &
              .AND. zs(igrd,jgrd,k) <= ht_defined ) THEN
          CALL spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd             &
                ,max_cld_snd,nz,iloc,jloc,k,cover,1.)
        END IF
      END DO
    ELSE
      DO k=2,k_ceil
        IF( wtcldcv(igrd,jgrd,k) /= 1.00                                &
              .AND. zs(igrd,jgrd,k) <= ht_defined ) THEN
          CALL spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd             &
                ,max_cld_snd,nz,iloc,jloc,k,cover,1.)
        END IF
      END DO
    END IF

    125     CONTINUE

  END DO       ! I=1,NOBSNG

  CALL mpiprocess_update(n_cld_snd,indexsng_wrk)
  CALL mpi_2dcr_collect(cld_snd,max_cld_snd,nz,n_cld_snd,indexsng_wrk)
  CALL mpi_2dcr_collect(wt_snd,max_cld_snd,nz,n_cld_snd,indexsng_wrk)

!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)                                                      &
    WRITE(6,*)' Num stations analyzed/cloud soundings = '               &
                            ,n_analyzed,n_cld_snd
!
!-----------------------------------------------------------------------
!
  istatus = 1
  RETURN
END SUBROUTINE insert_sao1
!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######        SUBROUTINE MODIFY_SOUNDING                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE modify_sounding (cf_modelfg,t_modelfg,hterain,               &
           t_sfc_k,i_in,j_in,k,nx,ny,nz,zs,                             &
           ht_base,ht_top,init,l_dry)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check if there is dry layer or inversion (If there is one,
!  ( the observed cloud layer will be cleared out).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  03/96   Based on the LAPS cloud analysis code, 1995
!
!  MODIFICATION HISTORY:
!
!  02/06/96  J. Zhang
!            Modified for ADAS grid. Added documents.
!
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!  INPUT:

  INTEGER :: init      ! flag value indicating if the layer is
                       ! within the cloud.

  INTEGER :: nx,ny               ! horizontal grid domain size
  INTEGER :: nz                  ! # of cloud analysis levels

  INTEGER :: i_in,j_in           ! i- and j- location in ADAS grid
  INTEGER :: k                   ! the cloud height level

  REAL :: cf_modelfg (nx,ny,nz)  ! first guess cloud cover field
  REAL :: t_modelfg (nx,ny,nz)   ! first guess temperature field
  REAL :: t_sfc_k (nx,ny)        ! surface temperature field
  REAL :: hterain(nx,ny)            ! height of the terrain
!
  REAL :: zs(nx,ny,nz)             ! cloud analysis heights
  REAL :: ht_base,ht_top
!
!-----------------------------------------------------------------------
!
!c    OUTPUT:

  LOGICAL :: l_dry                ! if there is inversion or dry layer?
!
!-----------------------------------------------------------------------
!
!c    LOCAL:

  INTEGER :: i,j
  REAL :: cf_model_base,t_model_base,t_subcloud
  REAL :: t_dry_adiabat,t_inversion_strength

  LOGICAL :: l_wait_for_base,l_cf,l_inversion

  SAVE l_wait_for_base,cf_model_base,t_model_base,                      &
       l_inversion,t_subcloud

!-----------------------------------------------------------------------
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  l_dry = .false.
  l_cf = .false.
!
!-----------------------------------------------------------------------
!
!    Find ARPS grid point nearest the SAO if it is out of bounds
!
!-----------------------------------------------------------------------
!
  i = MAX(MIN(i_in,nx-1),1)
  j = MAX(MIN(j_in,ny-1),1)

  IF (init == 1) THEN ! (below base )
                      ! reset to wait for the beginning of the layer
    l_wait_for_base = .true.
    l_inversion = .false.
    t_subcloud = t_modelfg(i,j,k)
    RETURN

  ELSE ! init = 0 (inside cloud layer)
       ! For k=1 and init=0, l_wait_base=.F., and there will be NO
       ! t_subcloud available.

    IF(l_wait_for_base .OR. k <= 2) THEN
!
!-----------------------------------------------------------------------
!
!      Set reference (just within cloud base)
!
!-----------------------------------------------------------------------
!
      l_wait_for_base = .false.
      cf_model_base = cf_modelfg(i,j,k)
      t_model_base = t_modelfg(i,j,k)

!          write(6,21) t_subcloud
!21        format( /' model T_subcloud = ',f7.2)
!            write(6,1)
!1           format(1x,21x,'  cf_fg','  t_fg ',' dt_inv'
!     :                 ,' lcf',' linv','  i   j   k ',' cldht')

    END IF
!
!-----------------------------------------------------------------------
!
!  Determine if the layer is dry or it has inversion.
!  (in either case, the cloud will be cleared out)
!
!-----------------------------------------------------------------------
!
    IF(.true.) THEN     ! Set inversion strength flag
!
      t_dry_adiabat = t_sfc_k(i,j)                                      &
                   -.0098 * (zs(i,j,k) - hterain(i,j))

      t_inversion_strength = t_modelfg(i,j,k) - t_dry_adiabat

      IF( ( (t_modelfg(i,j,k) > t_model_base)                           &
                            .OR.                                        &
           (k >= 2 .AND. t_modelfg(i,j,k) > t_subcloud) )               &
          .AND.                                                         &
           (t_modelfg(i,j,k) > 283.15)       & ! temp check
      .AND.                                                             &
           (t_inversion_strength > 4.)       & ! delta theta chk
      ) THEN
!
      l_inversion = .true.                  ! Inversion exists

!            write(6,2) cf_modelfg(i,j,k),t_modelfg(i,j,k)
!     :                 ,t_inversion_strength,l_cf,l_inversion
!     :                 ,i,j,k,nint(zs(i,j,k))
!
!2           format(' Inversion detected = ',3f7.2,2l4,3i4,i6)

    ELSE IF (cf_modelfg(i,j,k) < cf_model_base - 0.3   & ! cf search
        .AND. zs(i,j,k) - ht_base >= 500.) THEN
!
      l_cf = .true.           ! Dry layer exists

!                write(6,3) cf_modelfg(i,j,k),t_modelfg(i,j,k)
!     :               ,t_inversion_strength,l_cf,l_inversion
!     :               ,i,j,k,nint(zs(i,j,k))
!
!3               format(' Dry layer detected = ',3f7.2,2l4,3i4,i6)

    ELSE                          ! neither dry nor inversion
!
!-----------------------------------------------------------------------
!
!  If l_inversion and l_cf are not reset to .false. here, the WHOLE
!  cloud layer above a SINGLE level will be cleared out IF that one
!  level has inversion.
!  This is effective through the "save" statement at the beginning
!  of this subroutine.
!
!-----------------------------------------------------------------------
!
!            write(6,4) cf_modelfg(i,j,k),t_modelfg(i,j,k)
!     :              ,t_inversion_strength,l_cf,l_inversion
!     :              ,i,j,k,nint(zs(i,j,k))
!
!4           format(' model RH/T in cloud =',3f7.2,2l4,3i4,i6)

    END IF   ! inversion check

    IF( l_cf.OR.l_inversion ) THEN
      l_dry = .true.
    END IF

  END IF     ! .true. for dry-inversion check.

  END IF       ! INIT=1
!
!-----------------------------------------------------------------------
!
  RETURN
END SUBROUTINE modify_sounding

!
!##################################################################
!##################################################################
!######                                                      ######
!######              FUNCTION CLD_THK                        ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION cld_thk (ht_base,cover)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Obtain thickness of a cloud layer. The thickness is simply
!  set to 1 km for cloud below 7 km and 1.5km for those above
!  7 km.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  (Jian Zhang)
!  03/96
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    ht_base          ! the height of cloud base
!    cover            ! the cloud amount of the layer
!
!  OUTPUT:
!    cld_thk          ! the thickness of the cloud layer
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: ht_base,cover,cld_thk
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(ht_base > 7000.)THEN
    cld_thk = 1500.
  ELSE IF(ht_base > 1000.) THEN
    cld_thk = 1000.0+(ht_base-1000.)/12.
  ELSE
    cld_thk = MIN(1000.0,ht_base)
  END IF
  IF(cover > 0.01 .AND. cover < 0.5)THEN
    cld_thk = MIN(1000.,cld_thk)
  END IF

  RETURN
  END FUNCTION cld_thk

!
!##################################################################
!##################################################################
!######                                                      ######
!######              FUNCTION CLD_BASE                       ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION cld_base (k)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Obtain the base height of a cloud layer.
!  The base height is currently an artificial function of
!  the reported cloud layer index.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  (Jian Zhang)
!  04/17/96
!
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    k          ! the index of the reported cloud layer
!
!  OUTPUT:
!    cld_base   ! the base height of the cloud layer
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: cld_base
  INTEGER :: k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(k == 1) cld_base = 500.
  IF(k == 2) cld_base = 2000.
  IF(k == 3) cld_base = 6000.
  IF(k == 4) cld_base = 8000.
  IF(k == 5) cld_base = 10000.
  RETURN
  END FUNCTION cld_base





!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE INSERT_RADAR                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE insert_radar(nx,ny,nz,clddiag,hterain,zs,temp_3d,z_lcl       &
           ,ref_base1,ref_base2,hgt_thresh_ref                          &
           ,grid_ref,cldcv                                              &
           ,cloud_base,cloud_base_buf,l_unresolved                      &
           ,istatus)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Insert radar data into cloud grid
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  05/1996
!
!  MODIFICATION HISTORY:
!
!  05/08/96  J. Zhang
!            Modified for the ARPSDAS format. Added full
!            documentation.
!  07/26/96  J. Zhang
!            Added quality check to avoid out-bounds cloud.
!  03/27/97  J. Zhang
!            Updated for the official version of arps4.2.4
!  03/28/97  J. Zhang
!            Using the lifting condensation level as the cloud base
!            when there is no SAO cloud base.
!  09/01/97  J. Zhang
!            Document cleanup for ADASCLD version 25.
!            Using the lifting condensation level as the cloud base
!            when the radar echo top is below the analyzed SAO cloud
!            base
!  09/10/97  J. Zhang
!            Put lcl in the input argument list.
!  04/20/98  J. Zhang
!            Based on the arps4.3.3 version, abandon the cloud grid.
!            Using arps grid.
!  04/10/03  K. Brewster
!            Modify print statements to remove format numbers.
!            Renamed a few variables for consistency with ADAS.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!

  INTEGER :: nx,ny,nz
  !  INPUT/OUTPUT:
  REAL :: cldcv(nx,ny,nz)

  !  INPUT:

  INTEGER :: clddiag
  REAL :: temp_3d(nx,ny,nz)
  REAL :: zs(nx,ny,nz)
  REAL :: hterain(nx,ny)
  REAL :: z_lcl(nx,ny)      ! lifting condensation level

  REAL :: grid_ref(nx,ny,nz)

  REAL :: ref_base1         ! "significant" radar echo at lower levels
  REAL :: ref_base2         ! "significant" radar echo at upper levels
  REAL :: hgt_thresh_ref    ! height criteria for "significant" radar
                            ! echo thresholds
!  OUTPUT:
  INTEGER :: istatus
  REAL :: cloud_base(nx,ny)
  REAL :: cloud_base_buf(nx,ny) ! Lowest SAO/IR base within search radius
!
!  LOCAL:
  REAL :: radar_cover
  PARAMETER(radar_cover=1.0)
  REAL :: thresh_cvr   ! lower radar echo threshold for cloud filling
  PARAMETER (thresh_cvr = 0.2)


  LOGICAL :: l_below_base
  LOGICAL :: l_unresolved(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  LOGICAL :: cldprt
  REAL :: unlimited_ceiling,ref_thresh
  REAL :: zs_1d(nz)

  INTEGER :: i,j,k,kp1
  INTEGER :: icount_below,isearch_base,insert_count_tot
  INTEGER :: icount_radar_lvl,insert_count_lvl
  INTEGER :: i_out_bnd,i_lowecho

  INTEGER :: ips, ipe, jps, jpe

  LOGICAL :: verbose
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  verbose = .FALSE.

  IF (myproc == 0) WRITE(6,*) 'Inserting radar reflectivity '
  cldprt=(clddiag == 1)
!
!-----------------------------------------------------------------------
!
!  Calculate Cloud Bases
!
!-----------------------------------------------------------------------
!
  unlimited_ceiling = 200000.

  DO j = 1,ny-1
    DO i = 1,nx-1
      cloud_base(i,j) = unlimited_ceiling
!
      DO k = nz-1,1,-1
        IF(cldcv(i,j,k) < thresh_cvr .AND. cldcv(i,j,k+1) >= thresh_cvr) THEN
          cloud_base(i,j) = 0.5 * (zs(i,j,k) + zs(i,j,k+1))
        END IF
      END DO ! k
!
      IF (cloud_base(i,j) > z_lcl(i,j)) THEN
        cloud_base(i,j) = z_lcl(i,j)             ! using lcl
      END IF

      cloud_base_buf(i,j) = cloud_base(i,j)
      l_unresolved(i,j) = .false.

    END DO ! i
  END DO ! j

  icount_below = 0     ! tot # of rdr echo pts below cloud base
  isearch_base = 0     ! # of neighb cloud bases being succes. found
  insert_count_tot = 0 ! tot # of data pts inserted the radar_cover
!
!-----------------------------------------------------------------------
!
!  Essentially, this go downward to detect radar tops in time
!  to search for a new cloud base
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    ips = 2             ! Patch index
    ipe = nx-2
    jps = 2
    jpe = ny-2
    IF (loc_x == 1)       ips = 1         ! to cover the same domain as non-MPI runs below
    IF (loc_x == nproc_x) ipe = nx-1
    IF (loc_y == 1)       jps = 1
    IF (loc_y == nproc_y) jpe = ny-1
  ELSE
    ips = 1             ! Memory index for the whole domain
    ipe = nx-1
    jps = 1
    jpe = ny-1
  END IF

  i_out_bnd=0
  i_lowecho=0

  DO j = jps,jpe
    DO i = ips,ipe

      DO k=1,nz
        zs_1d(k) = zs(i,j,k)
      END DO

      icount_radar_lvl = 0 ! # of lvls with refl > ref_thresh
      insert_count_lvl = 0 ! # of cloud lvls inserted radar_cover

      DO k = nz-1,2,-1
        kp1 = MIN(k+1,nz)
!
!-----------------------------------------------------------------------
!
!  Define the "significant" reflectivity value based on height (to
!  eliminate ground clutters and/or other non-weather radar echoes).
!
!-----------------------------------------------------------------------
!
        IF ((zs(i,j,k)-hterain(i,j)) <= hgt_thresh_ref) THEN
          ref_thresh = ref_base1
        ELSE
          ref_thresh = ref_base2
        END IF

        IF(grid_ref(i,j,k) > ref_thresh) THEN
          icount_radar_lvl = icount_radar_lvl + 1

          l_below_base = .false.
!
!-----------------------------------------------------------------------
!
!  Test if we are at echo top
!
!-----------------------------------------------------------------------
!
          IF(k == nz-1 .OR. grid_ref(i,j,kp1) < ref_thresh) THEN
!
!-----------------------------------------------------------------------
!
!  Test if we are below the cloud base.
!
!-----------------------------------------------------------------------
!
            IF(zs(i,j,k) < cloud_base_buf(i,j)) THEN
!
!-----------------------------------------------------------------------
!
!  Radar Echo Top is below the analyzed cloud base.
!  Using the lifting condensation level as the modified cloud base
!  if it is lower than the current buffer value.
!
!-----------------------------------------------------------------------
!
              cloud_base_buf(i,j)=MIN(z_lcl(i,j),cloud_base_buf(i,j))

              IF(cloud_base_buf(i,j) < zs(i,j,k)) THEN
                isearch_base = isearch_base + 1
                IF(isearch_base < 50 .AND. verbose) THEN ! limit log output
                  WRITE(6,71) i,j,k,zs(i,j,k),cloud_base(i,j),cloud_base_buf(i,j)
                  71 FORMAT(' Rdr Top < Bse*gp=',3I3,' zs=',f7.0,       &
                            ' cld_b=',f7.0,'lcl=',f7.0,' Resolved')
                END IF

              ELSE ! Probably Unresolved base
                IF (verbose ) THEN
                  WRITE(6,72) i,j,k,zs(i,j,k),cloud_base(i,j),cloud_base_buf(i,j)
                  72 FORMAT(1X,' **** Prob Unresolved ****'/            &
                  &         1X,'Rdr Top < Bse*gp=',3I3,' zs=',f7.0,     &
                  &            ' cld_b=',f7.0,' lcl=',f7.0)
                END IF

                IF(cloud_base_buf(i,j) == unlimited_ceiling) THEN
                  l_unresolved(i,j) = .true.
                  WRITE(6,'(a,i4,i4,a)') 'Error, no LCL found for grid,',i,j, &
                          ' aborting from INSERTRAD...'
                  CALL arpsstop('No LCL found.',1)
                  !STOP
                END IF ! cloud_base(i,j).eq.unlimited_ceiling
                GO TO 600

              END IF ! cloud_base_buf(i,j) .lt. zs(i,j,k)

            END IF ! Blw Cld Base: zs.lt.cloud_base_buf(i,j)

          END IF ! At Echo Top: ref(k)>thr & (k.eq.nz .or. ref(kp1)<thr
!
!-----------------------------------------------------------------------
!
!  Loop through range of cloud grid levels for this model level
!
!-----------------------------------------------------------------------
!
          IF(zs(i,j,k) > cloud_base_buf(i,j)) THEN
                             ! Insert radar if we are above cloud base
            cldcv(i,j,k) = radar_cover
            insert_count_lvl = insert_count_lvl + 1
            insert_count_tot = insert_count_tot + 1
          ELSE ! Radar Echo below cloud base
            l_below_base = .true.
          END IF

          IF(l_below_base) THEN
            icount_below = icount_below + 1

            IF(icount_below <= 50 .AND. verbose) THEN
              WRITE(6,81) i,j,k,zs(i,j,k),cloud_base_buf(i,j)
              81 FORMAT(1X,'Rdr < Bse* g.p.:',3I3,                &
                           ' zs(i,j,k)=',f7.0,' cld_base=',f7.0)
            END IF

            GO TO 600
          END IF ! l_below_base =.true.

        ELSE ! grid_ref(i,j,k) <= ref_thresh
          IF(cldcv(i,j,k) > 0.4.AND.i_lowecho <= 25) THEN
            i_lowecho=i_lowecho+1
            IF(cldprt) WRITE(6,'(1x,a,f8.1,a,f6.2,a,3I3)')              &
            ' Reflect:',grid_ref(i,j,k),' <<< Cld cvr',cldcv(i,j,k),    &
            ' i,j,k=',i,j,k
          END IF
          GO TO 600
        END IF ! grid_ref(i,j,k) > ref_thresh ?

        IF (cldprt .AND. MOD(insert_count_tot,100) == 0 .AND. verbose) THEN
          WRITE(6,'(1x,a,i3,/,a,i5,a,i5,a,i5)')                         &
                'Inserted radar** k=',k,' nrdrl=',insert_count_lvl,     &
                ' n_instot=',insert_count_tot
        END IF
        600 CONTINUE

      END DO ! k
    END DO ! j
  END DO ! i

  CALL mpsumi(insert_count_tot,1)  ! it should be same as non-MPI runs
  IF (myproc == 0)                                                      &
    WRITE(6,*)' Total cloud grid points modified by radar = '           &
             ,insert_count_tot
  istatus=1

  RETURN
END SUBROUTINE insert_radar
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE INSERT_VIS                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE insert_vis (nx,ny,nz,zs,hterain,clouds_3d                    &
           ,albedo,cld_frac_vis_a                                       &
           ,vis_rad_thcvr,vis_rad_thdbz                                 &
           ,istat_radar,radar_ref_3d,ref_base,dbz_max_2d                &
           ,r_missing,sfc_sao_buffer,lvldbg,istatus)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To process satellite data for clouds and to modify the
!  3d cloud cover array.
!  Currently (10.7 micron) brightness temps are used with a
!  dummy call for the CO2 slicing method.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  04/26/96 (Based on LAPS code of 09/95)
!
!  MODIFICATION HISTORY:
!
!  04/26/96 (J. Zhang)
!           Modified for the ARPSDAS format.
!           Added full documentation.
!  09/02/97 (J. Zhang)
!           Added if statements to avoid putting in the reflectivity
!           threshold values in the data holes.
!  04/20/98 (J. Zhang)
!           Abandoned the cloud analysis grid, using the ARPS grid
!           instead.
!  04/10/03 K. Brewster
!           Modify print statements to remove format numbers.
!           Renamed a few variables for consistency with ADAS.
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
  INTEGER :: nx,ny,nz
  INTEGER, INTENT(IN) :: lvldbg
!  INPUT/OUTPUT:
  REAL :: clouds_3d(nx,ny,nz) ! cloud cover analysis
!
!  INPUT:
  REAL :: r_missing

  REAL :: zs(nx,ny,nz)
  REAL :: hterain(nx,ny)
!
!  INPUT:  parameters
  REAL :: vis_rad_thcvr          ! =0.2, defined in cloud_cv.f
  REAL :: vis_rad_thdbz          ! =10dbz, defined in cloud_cv.f
  REAL :: sfc_sao_buffer         ! =800m, defined in cloud_cv.f
!
!  INPUT: satellite visible channel data
  REAL :: albedo(nx,ny)
  REAL :: cld_frac_vis_a(nx,ny)
!
!  INPUT: radar data
  INTEGER :: istat_radar
  REAL :: dbz_max_2d(nx,ny)
  REAL :: radar_ref_3d(nx,ny,nz)
  REAL :: ref_base
!
!  OUTPUT:
  INTEGER :: istatus
!
!  LOCAL:
  INTEGER :: ih_alb(-10:20)
  INTEGER :: ih_cf_sat(-10:20)           ! Histogram for VIS cf array
  INTEGER :: ih_cfin(-10:20)             ! Histogram for input cf array
  INTEGER :: ih_cfout(-10:20)            ! Histogram for output cf array
  INTEGER :: ih_cfin_out(-10:20,-10:20)  ! Histogram for the comparison
                                       ! between input cloud cover
                                       ! array and modified cf array.
  INTEGER :: ih_cfin_sat(-10:20,-10:20)  ! Histogram for the comparison
                                       ! between input cloud cover
                                       ! array and sat. VIS cf array.
  INTEGER :: ih_cmaxin_sat(-10:20,-10:20)
  INTEGER :: ih_cmaxout_sat(-10:20,-10:20)

!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL :: viscorfrc
  PARAMETER (viscorfrc=0.29)

  INTEGER :: i,j,k,kl
  INTEGER :: n_missing_albedo,n_vis_mod
  INTEGER :: i_sat,i_in,i_out,i_cmaxin,i_cmaxout
  REAL    :: diffin,diffout,diffin_sum,diffout_sum,                     &
             diffin_sumsq,diffout_sumsq
  REAL    :: cld_frac_vis,cld_frac_in,cld_frac_out
  REAL    :: cushion
  REAL    :: cmaxin,cmaxout
  INTEGER :: iblank_radar,iset_vis
  REAL    :: r_present
  LOGICAL :: l_prt

  INTEGER :: nylg        ! Add by Y. Wang for MPI diagnostic outputs
  INTEGER :: jmid,jmidproc
  LOGICAL :: verbose = .FALSE.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus=0

  IF (lvldbg > 90) verbose = .TRUE.

  nylg = (ny-3)*nproc_y + 3
  jmid = nylg/2
  jmidproc = (jmid-2) / (ny-3) + 1
  IF (loc_y == jmidproc) THEN
    jmid = MOD((jmid-2),(ny-3)) + 2   ! Local index for global middle point
  ELSE
    jmid = -999   ! a unreasonable value
  END IF
!
!-----------------------------------------------------------------------
!
!  Initialize all histogram arrays.
!
!-----------------------------------------------------------------------
!
  DO i=-10,20
    ih_alb(i) = 0
    ih_cf_sat(i) = 0
    ih_cfin(i) = 0
    ih_cfout(i) = 0
    DO j=-10,20
      ih_cfin_out(i,j) = 0
      ih_cfin_sat(i,j) = 0
      ih_cmaxin_sat(i,j) = 0
      ih_cmaxout_sat(i,j) = 0
    END DO
  END DO
!
  n_missing_albedo = 0
  n_vis_mod = 0        ! # of data points modified by VIS data

  diffin_sum  = 0.     ! sum of diff. between column max. INPUT cld
                        ! fraction and the satellite vis cld frac'n.
  diffout_sum = 0.     ! sum of diff. between column max. OUTPUT cld
                        ! fraction and the satellite vis cld frac'n.
  diffin_sumsq  = 0.
  diffout_sumsq = 0.
!
!-----------------------------------------------------------------------
!
!  Horizontal array loop
!
!-----------------------------------------------------------------------
!
  l_prt = .false.
  DO i = 1,nx-1
    DO j = 1,ny-1
!
      kl  = nint(albedo(i,j)*10.)
      kl  = MIN(MAX(kl,-10),20)
      ih_alb(kl) = ih_alb(kl) + 1

      IF(cld_frac_vis_a(i,j) /= r_missing) THEN

        cld_frac_vis = cld_frac_vis_a(i,j)

        i_sat = nint(cld_frac_vis*10.)
        i_sat = MIN(MAX(i_sat,-10),20)
        ih_cf_sat(i_sat) = ih_cf_sat(i_sat) + 1  ! # of data points
                                  ! in a VIS cf value catagory
!
!-----------------------------------------------------------------------
!
!  Make sure satellite cloud fraction is between 0 and 1
!
!-----------------------------------------------------------------------
!
        IF(cld_frac_vis < 0.0) cld_frac_vis = 0.0
        IF(cld_frac_vis > 1.0) cld_frac_vis = 1.0

        i_sat = nint(cld_frac_vis*10.)

        cmaxin = 0.
        cmaxout = 0.

        iblank_radar = 0
        iset_vis = 0

        DO k = 2,nz-1

!          l_out(1)=k.eq.19.and.i.le.17.and.i.ge.2
!     :                    .and.j.le.43.and.j.ge.36
!          l_out(3)=k.eq.26.and.i.le.10.and.i.ge.6
!     :                    .and.j.le.26.and.j.ge.30
!          l_prt=l_out(1).or.l_out(3)

          cld_frac_in = MIN(clouds_3d(i,j,k),1.0) ! save input cf value
!
!-----------------------------------------------------------------------
!
!  Modify the cloud field with the vis input - allow .3 vis err?
!  The scheme:
!  Above the sao buffer layer, the previoud cloud analysis is
!  retained if it is within cld_frac_vis + 0.3.
!  Below the layer, the cloud analysis must be smaller than
!  cld_frac_vis. If not, reset it to cld_frac_vis.
!
!-----------------------------------------------------------------------
!
          IF(zs(i,j,k) > hterain(i,j)+sfc_sao_buffer) THEN
            cushion = 0.3     ! Allow 0.3 error in VIS cld cvr.
          ELSE
            cushion = 0.0     ! VIS cld cover is the max. allowed.
          END IF

          IF(cld_frac_vis < viscorfrc .AND.                             &
                ((clouds_3d(i,j,k)-cld_frac_vis) > cushion)) THEN
            cld_frac_out = cld_frac_vis     ! VIS cf value
!
            IF(l_prt) THEN
              WRITE(6,621) i,j,k,clouds_3d(i,j,k)                       &
                           ,radar_ref_3d(i,j,k),dbz_max_2d(i,j)         &
                           ,cld_frac_vis_a(i,j)                         &
                           ,cld_frac_out
              621  FORMAT(1X,3I3,' oldcld=',f5.2,' ref=',f8.1,' dbz_m=',f8.1 &
                           ,' viscld=',f5.2,' newcld=',f5.2)
            END IF

!
!-----------------------------------------------------------------------
!
!  Determine if we need to reconcile VIS with radar
!
!-----------------------------------------------------------------------
!
            IF(istat_radar == 1 .AND. dbz_max_2d(i,j) /= r_missing      &
                  .AND. dbz_max_2d(i,j) > ref_base) THEN
                                                     ! Valid radar echo

              IF(cld_frac_out < vis_rad_thcvr) THEN

                IF(dbz_max_2d(i,j) < vis_rad_thdbz) THEN
                                  ! Blank out Radar, Normal VIS Clearing
                  iblank_radar = 1
                ELSE            ! dbz_max_2d(i,j) >= vis_rad_thdbz
                  cld_frac_out = vis_rad_thcvr     ! use radar cf
                  iset_vis = 1
                END IF

              END IF
            END IF

            IF(cld_frac_in - cld_frac_out > .01) THEN
              n_vis_mod = n_vis_mod + 1
            END IF

            clouds_3d(i,j,k) = cld_frac_out   ! Modify the output
          ELSE  ! clouds_3d - cld_frac_vis .le. cushion
            cld_frac_out = cld_frac_in        ! previous cloud analysis
          END IF  ! clouds_3d - cld_frac_vis .gt. cushion
!
!-----------------------------------------------------------------------
!
!  Update Histograms
!
!  Check to make sure we don't to outside the bounds...every once in a
!  while it happens.  Model after what is done with "i_sat".
!
!-----------------------------------------------------------------------
!
          i_in = nint(cld_frac_in*10.)
          i_in = MIN(MAX(i_in,-10),20)
          ih_cfin(i_in) = ih_cfin(i_in) + 1

          i_out = nint(cld_frac_out*10.)
          i_out = MIN(MAX(i_out,-10),20)
          ih_cfout(i_out) = ih_cfout(i_out) + 1

          ih_cfin_sat(i_in,i_sat)                                       &
                         = ih_cfin_sat(i_in,i_sat) + 1

          ih_cfin_out(i_in,i_out)                                       &
                         = ih_cfin_out(i_in,i_out) + 1

          cmaxin  = MAX(cmaxin,cld_frac_in)   !col.max of inp cldcvr
          cmaxout = MAX(cmaxout,cld_frac_out) !col.max of outp cldcvr

        END DO ! k
!
!-----------------------------------------------------------------------
!
!  Reconcile VIS with radar
!
!-----------------------------------------------------------------------
!
        IF(iblank_radar == 1) THEN ! NO VIS / WEAK ECHO
!
!-----------------------------------------------------------------------
!
!  Blank out radar column for this grid point
!
!-----------------------------------------------------------------------
!
!         WRITE(6,'(a,2i3,f8.2,f8.1,f8.2)')                             &
!           ' Vis_Rdr - Blank out radar: cvr/refl/vis << ',             &
!           i,j,cmaxout,dbz_max_2d(i,j),cld_frac_vis
          IF (dbz_max_2d(i,j) > ref_base) dbz_max_2d(i,j) = ref_base
          DO kl = 1,nz
            radar_ref_3d(i,j,kl)=min(radar_ref_3d(i,j,kl),ref_base)
          END DO ! kl

        ELSE IF (iset_vis == 1) THEN ! NO VIS / STRONG ECHO
!
!-----------------------------------------------------------------------
!
!  Cloud cvr has been reset to threshold value above VIS
!
!-----------------------------------------------------------------------
!
!         WRITE(6,'(a,2i4,f5.2,f8.1,2f5.2)') &
!             ' VIS_RDR - Reset vis: cvr/dbz/vis/thr << ', &
!             i,j,cmaxout,dbz_max_2d(i,j),cld_frac_vis,vis_rad_thcvr

        END IF ! iblank_radar=1

        IF(j == jmid .AND. (cmaxin-cmaxout) > 0.1) THEN
          WRITE(6,'(1x,a,2i3,a,f5.2,a,f5.2)') &
                'Vismod',i,j,' colcfi=',cmaxin,' colcfo=',cmaxout
        END IF

        i_cmaxin  = nint(cmaxin*10.)
        i_cmaxout = nint(cmaxout*10.)
        ih_cmaxin_sat(i_cmaxin,i_sat)                                   &
                          = ih_cmaxin_sat(i_cmaxin,i_sat) + 1
        ih_cmaxout_sat(i_cmaxout,i_sat)                                 &
                          = ih_cmaxout_sat(i_cmaxout,i_sat) + 1

        diffin  = cmaxin  - cld_frac_vis
        diffout = cmaxout - cld_frac_vis
        diffin_sum  = diffin_sum  + diffin
        diffout_sum = diffout_sum + diffout
        diffin_sumsq  = diffin_sumsq  + diffin**2
        diffout_sumsq = diffout_sumsq + diffout**2

      ELSE     ! cld_frac_vis_a(i,j) = r_missing
        n_missing_albedo =  n_missing_albedo + 1

      END IF    ! cld_frac_vis_a(i,j) .ne. r_missing

    END DO ! i
  END DO ! j

  IF (Verbose) THEN

    WRITE(6,'(/a,i10)') ' N_MISSING_ALBEDO =',n_missing_albedo
    WRITE(6,'(a,i10)') ' N_VIS_MOD =',n_vis_mod
    WRITE(6,'(/a)')'           HISTOGRAMS'
    WRITE(6,'(a)')'   I    Alb  CFsat   CFin  CFout'
    DO i = -5,15
      WRITE(6,'(1X,i3,4i7)')                                              &
            i,ih_alb(i),ih_cf_sat(i),ih_cfin(i),ih_cfout(i)
    END DO

    WRITE(6,'(/a/)') ' Input vs. Satellite Cloud Fraction Histogram'
    WRITE(6,'(a)') ' X-axis: input cld frac, Y-axis: Satellite cld frac'
    WRITE(6,'(4x,11i6)') (i,i=0,10)
    DO j = 0,10
      WRITE(6,'(1x,i3,11i6)') j,(ih_cfin_sat(i,j),i=0,10)
    END DO

    WRITE(6,'(a)')'  Input vs. Output Cloud Fraction Histogram'
    WRITE(6,'(a)')' X-axis: input cld frac, Y-axis: output cld fraction'
    WRITE(6,'(4x,11i6)') (i,i=0,10)
    DO j = 0,10
      WRITE(6,'(1x,i3,11i6)') j,(ih_cfin_out(i,j),i=0,10)
    END DO

    WRITE(6,'(/a)') ' Column Max Input vs. Satellite Fraction Histogram'
    WRITE(6,'(/a)') ' X-axis: column max.input cf, Y-axis: col.max.sat.cf.'
    WRITE(6,'(4x,11i6)') (i,i=0,10)
    DO j = 0,10
      WRITE(6,'(1x,i3,11i6)') j,(ih_cmaxin_sat(i,j),i=0,10)
    END DO

    WRITE(6,'(/a)')'  Column Max Output vs. Satellite Fraction Histogram'
    WRITE(6,'(/a)')' X-axis: column max.output cf, Y-axis: col.max.sat.cf.'
    WRITE(6,'(4x,11i6)') (i,i=0,10)
    DO j = 0,10
      WRITE(6,'(1x,i3,11i6)') j,(ih_cmaxout_sat(i,j),i=0,10)
    END DO

    r_present = nx*ny - n_missing_albedo
    IF(r_present > 0.) THEN ! write stats
      WRITE(6,'(a,2f8.3)') ' VIS STATS: Mean/RMS input residual  = ', &
         diffin_sum/r_present,SQRT(diffin_sumsq/r_present)
      WRITE(6,'(a,2f8.3)')  ' VIS STATS: Mean/RMS output residual = ', &
         diffout_sum/r_present,SQRT(diffout_sumsq/r_present)
    END IF
  END IF

  istatus = 1

  RETURN
END SUBROUTINE insert_vis
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE INSERT_IR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE insert_ir (nx,ny,nz,rlat,rlon,rland_frac,cvr_snow,           &
                      hterain,zs,z_lcl,temp_3d,p_3d,                    &
                      t_sfc_k,t_gnd_k,cldcv_sao,                        &
                      solar_alt,solar_ha,solar_dec,                     &
                      isatir,t11mu,cldtop_m_ir,cldtop_m,                &
                      dx,dy,sfc_sao_buffer,                             &
                      cldcv,tem1,lvldbg,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To process satellite data for clouds and to modify the
!  3d cloud cover array.
!  Currently 10,7 micron brightness temps are used.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  (Jian Zhang)
!  04/1996  Based on the LAPS cloud analysis code (Steve Albers,
!           09/1995)
!
!  MODIFICATION HISTORY:
!
!  04/26/96  J. Zhang
!            Modified for the ARPSDAS format.
!            Added full documentation.
!  07/19/96  J. Zhang
!            Added the quality check part to deal with the tall
!            and steep cloud towers.
!  07/24/96  J. Zhang
!            Added the quality check for cases when the cloud top
!            is higher than the model top boundary
!  11/01/96  J. Zhang
!            Added the rawinsonde data for determing the cloud
!            top height for the low clouds.
!  03/18/97  J. Zhang
!            Cleanup the code and implemented for the official
!            arps4.2.4 version
!  08/06/97  J. Zhang
!            Change adascld24.inc to adascld25.inc
!  09/09/97  J. Zhang
!            Fixed a bug when calling COMPARE_RAD. Further cleanup.
!  09/10/97  J. Zhang
!            Lifting condensation level is now an input argument.
!            Change adascld25.inc to adascld26.inc
!  04/20/98  J. Zhang
!            Abandoned the cloud analysis grid, using the ARPS grid
!            instead.
!  2002-06-12 G. Bassett
!            Redid the way t11mu_cold & ht_sao_base are computed so that no
!            grid points in the area concerned are skipped.
!  03/11/03  Keith Brewster, CAPS
!            Modified some code to streamline and new ir calibration.
!            Also changed some variable names to be consistent with the
!            current GOES imagers.  Changed missing value flags to ease
!            printing and changed some print formats.
!  03/31/07  Keith Brewster, CAPS
!            Moved ir_cold calculation to mci2arps and is not read
!            from the ir file in cloud_cv and passed into this subroutine
!            Changed array t10mu to 3-d to store both IR and cold-filtered
!            IR in levels 1 and 2, respectively.
!
!  04/12/07  Yunheng Wang, CAPS
!            Fixed for searching nearby SAO cld layers in MPI mode.
!
!  09/13/07  Kevin W. Thomas CAPS
!            In the non-MPI case, "cldcv_sao" and "cldcv" share the same
!            memory, which is very bad.  "cldcv_sao" know always goes to
!            a temporary array.
!            Fix a subscript problem that could cause a wrong "zs" value
!            to be assigned to a point.
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

  INCLUDE 'phycst.inc'
  INCLUDE 'adas.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz            ! ARPS grid size
  INTEGER, INTENT(IN) :: lvldbg
!  INCLUDE: ( from adas.inc)
!  real r_missing
!
!  INPUT/OUTPUT:
  REAL :: cldcv(nx,ny,nz)        ! 3D cldcv analx after useing sat.data
  REAL :: tem1(nx,ny,nz)         ! work variable
!
!  INPUT:   general
  REAL :: dx,dy
!
!  INPUT:  for model grid geometry
  REAL :: hterain(nx,ny)
  REAL :: rlat(nx,ny),rlon(nx,ny) ! Lat. and lon. of the grid pt.
  REAL :: zs(nx,ny,nz)            ! Phys. ht (m) at model grid point.
  REAL :: z_lcl(nx,ny)            ! Lifting condensation level (MSL)
  REAL :: temp_3d(nx,ny,nz)       ! temp. analysis on model grid
  REAL :: p_3d(nx,ny,nz)          ! pres. analysis on model grid
!
!  INPUT: for model grid sfc characteristic fields
  REAL :: rland_frac(nx,ny)  ! land coverage fraction at a g.p.
  REAL :: cvr_snow(nx,ny)         ! snow coverage at a grid point.
  REAL :: solar_alt(nx,ny)        ! solar altitude angle
  REAL :: solar_ha(nx,ny)         ! solar hour angle
  REAL :: solar_dec
!
!  INPUT: for model grid analysis fields
  REAL :: t_sfc_k(nx,ny)          ! sfc air temp.
!
!  INPUT:  for SAO
  REAL :: sfc_sao_buffer          ! No clearing cloud by sat. below
                                  ! this ht.(m AGL) (hence letting
                                  ! SAOs dominate)
  REAL, TARGET, INTENT(IN) :: cldcv_sao(nx,ny,nz) ! 3D Cloud cover array
!
!  OUTPUT:
  INTEGER :: istatus
  INTEGER :: isatir(nx,ny)
  REAL :: t_gnd_k(nx,ny)          ! ground skin temp.
  REAL :: t11mu(nx,ny,2)         ! Satell 10.7 mu IR bright. temp.
                                 ! Index 1: raw IR, 2: cold-filtered
  REAL :: cldtop_m(nx,ny)         ! Adj. Sat. cloud top height (m)
  REAL :: cldtop_m_ir(nx,ny)     ! Sat. cld top ht (m) from 10.7 mu IR
!
!  Local: factors and parameters
  REAL :: sfc_ir_buffer           !No adding cld by sat. below this ht.
                                  ! (m AGL) (hence letting SAOs dominate)
  PARAMETER (sfc_ir_buffer = 3000.)

  REAL :: thk_def                 !Def. cld thkness inserted by sat.data
  PARAMETER (thk_def = 1500.)

  REAL :: thr_sao_cvr             ! cldcvr thres. in evaluating presence
                                  ! of SAO cld layers
  PARAMETER (thr_sao_cvr = 0.1)

  REAL :: thresh2                 ! Thresh. for IR cloud detection
  PARAMETER (thresh2 = 3.5)       ! (Tsfc - T_IR) Were 5K and 15K also
!
!  LOCAL: for 10.7 micron IR temp.
  REAL :: t11mu_est               ! estimated 10.7 mu IR temp. by using
                                  ! cld cvr analysis
!  LOCAL: for interpolations
  REAL :: zs_1d(nz)
  REAL :: t_1d(nz)  ! 3D temperature analysis on model grid
  REAL :: p_1d(nz)  ! pres. analysis on model grid
!
!  LOCAL: for satellite cloud presence and cloud top height
  LOGICAL :: l_t11mu,l_cloud_present
  REAL :: cldtop_m_old
  INTEGER :: nlyr               ! for iterative adj. of cld cvr field
!
!  FUNCTIONS:
  REAL :: t_ground_k
  REAL :: temp_to_rad
  REAL :: rad_to_temp
  REAL :: ir_cover
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,i_sao,j_sao,isat
  INTEGER :: il,ih,jl,jh,ii,jj
  INTEGER :: n_no_sao1,n_no_sao2,n_no_sao3
  INTEGER :: i_reset_base

  REAL :: htbase,htbase_init,ht_cloud
  PARAMETER (ht_cloud = 2000.0)

  INTEGER :: idelt,jdelt,idelt_max,jdelt_max,idelt_index,jdelt_index

  REAL :: sat_cover,temp_thresh,cldcv_above,cover,cover_new,buffer
  REAL :: tmin

  LOGICAL :: verbose = .FALSE.

  REAL              :: ht_sao_base_min
  REAL, ALLOCATABLE :: ht_sao_base (:,:)
  REAL, ALLOCATABLE :: cldcv_sao_wrk(:,:,:)
  REAL, ALLOCATABLE :: zs_wrk(:,:,:)

  INTEGER           :: isrcbgn, isrcend, jsrcbgn, jsrcend
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (lvldbg > 90) verbose = .TRUE.

!-----------------------------------------------------------------------
!
!  Print some 10.7 micron IR brightness temp. samples.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,'(a/)') 'Insert satellite IR data routines called.'

  il=MAX(1,nx/2-5)
  ih=MIN(nx,nx/2+5)
  jl=MAX(1,ny/2-5)
  jh=MIN(ny,ny/2+5)

  IF ( verbose ) THEN    ! NOTE: only local diagnostic in processor 0
    WRITE(6,'(a)') ' 10.7-micron IR brightness temp values in the mid-domain:'
    WRITE(6,'(/4x,11(2x,i2,2x))') (i,i=il,ih)
    WRITE(6,'(1x,i3,11f6.1)') (j,(t11mu(i,j,1),i=il,ih),j=jh,jl,-1)
  END IF

!
!-----------------------------------------------------------------------
!
!  Define a box for searching SAO cloud base
!  ht_sao_base flags: -99999. => not computed yet, 99999. => no sao base.
!
!-----------------------------------------------------------------------
!
  idelt_max = nint(50000. / dx)
  jdelt_max = nint(50000. / dy)

  idelt_max = MIN(nx-3,idelt_max)  ! We only want to change message with
  jdelt_max = MIN(ny-3,jdelt_max)  ! neighbors
!
!-----------------------------------------------------------------------
!
!  Array "cldcv_sao" is an input array with SAO info.  "cldcv" is an
!  output array that initially contains the contents of "cldcv_sao".
!
!  Based on the code in this file, these are two different arrays.
!  Looking at the caller in "cloud_cv.f90", the user will notice that
!  variable "clouds_3d" is used in both places, meaning that in the
!  subroutine, "cldcv_sao" and "cldcv" both reference the same location
!  in memory.  This isn't good.  Too fix this, we put the "cldcv_sao"
!  data in a temporary array, then no longer user "cldcv_sao".
!
!-----------------------------------------------------------------------
!

  IF (mp_opt > 0) THEN
    isrcbgn = 1-idelt_max     ! Extended arrays to be used for searching SAO cloud base
    isrcend = nx+idelt_max
    jsrcbgn = 1-jdelt_max
    jsrcend = ny+jdelt_max

    ALLOCATE(cldcv_sao_wrk(isrcbgn:isrcend,jsrcbgn:jsrcend,nz), STAT = istatus)
    CALL check_alloc_status(istatus, "insert_ir:cldcv_sao_wrk")
    ALLOCATE(zs_wrk(isrcbgn:isrcend,jsrcbgn:jsrcend,nz), STAT = istatus)
    CALL check_alloc_status(istatus, "insert_ir:zs_wrk")

    CALL mpfillextarray3d(cldcv_sao,nx,ny,nz,idelt_max,jdelt_max,       &
               cldcv_sao_wrk,ebc,wbc,sbc,nbc,0,istatus)
    CALL mpfillextarray3d(zs,       nx,ny,nz,idelt_max,jdelt_max,       &
               zs_wrk,       ebc,wbc,sbc,nbc,0,istatus)

    IF (loc_x == 1) THEN         ! west boundary
      isrcbgn = 1
    END IF

    IF (loc_x == nproc_x) THEN   ! east boundary
      isrcend = nx
    END IF

    IF (loc_y == 1) THEN         ! south boundary
      jsrcbgn = 1
    END IF

    IF (loc_y == nproc_y) THEN   ! north boundary
      jsrcend = ny
    END IF

  ELSE
    isrcbgn = 1
    isrcend = nx
    jsrcbgn = 1
    jsrcend = ny

    ALLOCATE(cldcv_sao_wrk(isrcbgn:isrcend,jsrcbgn:jsrcend,nz), STAT = istatus)
    CALL check_alloc_status(istatus, "insert_ir:cldcv_sao_wrk")
    ALLOCATE(zs_wrk       (isrcbgn:isrcend,jsrcbgn:jsrcend,nz), STAT = istatus)
    CALL check_alloc_status(istatus, "insert_ir:zs_wrk")
    cldcv_sao_wrk(:,:,:) = cldcv_sao(:,:,:)
    zs_wrk       (:,:,:) = zs       (:,:,:)
  END IF

  ALLOCATE(ht_sao_base  (isrcbgn:isrcend,jsrcbgn:jsrcend),    STAT = istatus)
  CALL check_alloc_status(istatus, "insert_ir:ht_sao_base")

  ht_sao_base(:,:) = -99999.

  isrcend = isrcend -1   ! because ht_sao_base is a mass variable
  jsrcend = jsrcend -1

!
!-----------------------------------------------------------------------
!
!  First guess conversion from cloud height grid to ADAS grid
!  This has to err slightly on the high side
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  QC check on 10.7 micron brightness temps
!
!-----------------------------------------------------------------------
!
  DO j = 1,ny
    DO i = 1,nx
      IF(t11mu(i,j,1) < 173. .OR. t11mu(i,j,1) > 350.) t11mu(i,j,1)=-999.
      IF(t11mu(i,j,2) < 173. .OR. t11mu(i,j,2) > 350.) t11mu(i,j,2)=t11mu(i,j,1)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate ground temperature
!  This is air temperature over water and a radiation balanced condx
!  over land.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,*) 'Computing ground skin temperature'
  DO j = 1,ny-1
    DO i = 1,nx-1
      IF ( rland_frac(i,j) >= 0.5) THEN
        t_gnd_k(i,j) = t_ground_k(t_sfc_k(i,j),solar_alt(i,j)           &
                                  ,solar_ha(i,j),solar_dec              &
                                  ,rlat(i,j),cvr_snow(i,j)              &
                                  ,r_missing,i,j,nx,ny)
      ELSE ! water environment
        t_gnd_k(i,j) = t_sfc_k(i,j)
      END IF
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Initializes routine
!
!-----------------------------------------------------------------------
!
  n_no_sao1 = 0     ! # of g.p with sat.cld but no SAO cld at the grid pt
  n_no_sao2 = 0     ! # of g.p with sat.cld but no SAO cld in the box
  n_no_sao3 = 0     ! # of pts with sat.cldtop lower than SAO cld base.
  i_reset_base = 1  ! # of pts that sate. ceiling being resetted
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,*) 'Find cloud_top for each grid point.'

  DO j=1,ny-1
    DO i=1,nx-1
      IF(t11mu(i,j,1) > 0.) THEN
        isat=isatir(i,j)
        DO k=1,nz
          zs_1d(k) = zs(i,j,k)
          t_1d(k)  = temp_3d(i,j,k)
          p_1d(k)  = p_3d(i,j,k)
        END DO
        tmin=temp_3d(i,j,2)
        DO k=2,nz-1
          tmin=MIN(tmin,temp_3d(i,j,k))
        END DO
!
!-----------------------------------------------------------------------
!
!  Calculate cloud top height from 10.7-micron method
!
!-----------------------------------------------------------------------
!
        CALL cloud_top(nx,ny,nz,                                        &
                    i,j,zs_1d,t_1d,p_1d,z_ref_lcl,z_lcl(i,j),           &
                    hterain(i,j),tmin,r_missing,                        &
                    t11mu(i,j,1),t_gnd_k(i,j),thresh2,                  &
                    cldtop_m_ir(i,j),cldtop_m(i,j),                     &
                    l_t11mu,l_cloud_present,sat_cover)
!
!-----------------------------------------------------------------------
!
!  Modify those levels where the satellite shows warmer than the
!  calculated brightness temp from the analysis of SAO/Pireps
!
!-----------------------------------------------------------------------
!
        DO k = nz-1,1,-1

          IF (cldcv(i,j,k) > 0.04) THEN ! Efficiency test
!
!-----------------------------------------------------------------------
!
!  Estimate the brightness temperature using the cloud cover
!  field.
!
!-----------------------------------------------------------------------
!
            IF (cldcv(i,j,k) /= r_missing) THEN
              t11mu_est = temp_to_rad(isat,temp_3d(i,j,k))                &
                  * cldcv(i,j,k) + temp_to_rad(isat,t_gnd_k(i,j))       &
                  * (1.-cldcv(i,j,k))
              t11mu_est = rad_to_temp(isat,t11mu_est)  ! estimated IR temp.

!                t11mu_est_a = temp_3d(i,j,k) * cldcv(i,j,k)
!     :             + t_gnd_k(i,j) * (1.-cldcv(i,j,k))

            ELSE
              t11mu_est = t_gnd_k(i,j)
            END IF
!
!-----------------------------------------------------------------------
!
!  Test if clouds detected by SAO/PIREP should have been
!  detected by the satellite (if satellite is warmer than analysis)
!
!-----------------------------------------------------------------------
!
            IF (t11mu(i,j,1) > t11mu_est) THEN !seems too much SAO cld

              !Don't touch points within buffer of surface
              IF (zs(i,j,k)-hterain(i,j) > sfc_sao_buffer) THEN

                !Does satellite still imply at least some cloud?
                IF(t_gnd_k(i,j)-t11mu(i,j,1) > thresh2) THEN ! Some cloud
                  IF(cldcv(i,j,k) > 0.9) THEN
                    cldcv(i,j,k)=.01           ! Lower top of solid cld
                  ELSE                         ! Cover < 0.9, correct it
                    cldcv(i,j,k) = ir_cover(isat,t11mu(i,j,1)         &
                        ,t_gnd_k(i,j),temp_3d(i,j,k))

                    cldcv(i,j,k) = MAX(0.0,MIN(1.0,cldcv(i,j,k)))
                  END IF

                ELSE ! IR temp nearly matches grnd, clear it
!
!-----------------------------------------------------------------------
!
!           Insure that "black (or grey) stratus" is not present
!
!-----------------------------------------------------------------------
!
                  temp_thresh = MIN(t_gnd_k(i,j),t_sfc_k(i,j)-10.0)
                  IF(temp_3d(i,j,k) < temp_thresh)THEN
                    cldcv(i,j,k)=.01 ! not in inversion, clear it out
                  END IF

                END IF ! IR signature present
              END IF ! cldht-hterain.gt.sfc_sao_buffer is .true.
            END IF ! t11mu(i,j).gt.t11mu_est
          END IF ! Current Cloud Cover is significant (> .04)
        END DO ! k (for clearing clouds)
!
!-----------------------------------------------------------------------
!
!  Insert satellite clouds
!
!-----------------------------------------------------------------------
!
        IF(l_cloud_present) THEN ! Insert satellite clouds
!
!-----------------------------------------------------------------------
!
!  Set initial satellite cloud base.
!
!-----------------------------------------------------------------------
!
          htbase_init=cldtop_m(i,j) - thk_def
          htbase = htbase_init
!
!-----------------------------------------------------------------------
!
!  Locate lowest SAO cloud base
!
!-----------------------------------------------------------------------
!
         IF (ht_sao_base(i,j) < -90000.) THEN ! haven't tried yet

            ht_sao_base(i,j) = 99999.
            cldcv_above = cldcv_sao_wrk(i,j,nz-1)

            DO k=nz-2,1,-1

              IF (cldcv_sao_wrk(i,j,k) <= thr_sao_cvr                   &
                    .AND. cldcv_above > thr_sao_cvr) THEN
                ht_sao_base(i,j) = zs(i,j,k+1)
              END IF

              cldcv_above = cldcv_sao_wrk(i,j,k)
            END DO ! k

          ENDIF

          i_sao = i
          j_sao = j
          ht_sao_base_min = ht_sao_base(i,j)

          IF (ht_sao_base(i,j) > 90000.) THEN ! Srch for nearby SAO cld lyrs
                                       ! because the sat. says cld
                                       ! and the SAO doesn't
            n_no_sao1 = n_no_sao1 + 1 ! Sat.cld but no SAO cld at same pt.

            ht_sao_base_min = -99999.
            idelt_index = 1
            DO ! search for neighboring ht_sao_base values
              IF (idelt_index > idelt_max .OR. ht_sao_base_min > -90000.) EXIT

              DO jj = j-idelt_index,j+idelt_index
                DO ii = i-idelt_index,i+idelt_index,2*idelt_index
!                  IF(ii >= 1.AND.ii <= nx-1 .AND. jj >= 1.AND.jj <= ny-1 ) THEN
                  IF(ii >= isrcbgn .AND. ii <= isrcend .AND.            &
                     jj >= jsrcbgn .AND. jj <= jsrcend ) THEN   ! bounds check
                    IF (ht_sao_base(ii,jj) < -90000. ) THEN
                      ! ht_sao_base(ii,jj) not computed yet
                      ht_sao_base(ii,jj) = 99999.
                      cldcv_above = cldcv_sao_wrk(ii,jj,nz-1)
                      DO k = nz-2,1,-1
                        IF (cldcv_sao_wrk(ii,jj,k) <= thr_sao_cvr       &
                              .AND. cldcv_above > thr_sao_cvr) THEN
!                         ht_sao_base(ii,jj) = zs(i,j,k+1)
                          ht_sao_base(ii,jj) = zs_wrk(ii,jj,k+1)
                        END IF
                        cldcv_above = cldcv_sao_wrk(ii,jj,k)
                      END DO ! k
                    END IF ! compute ht_sao_base(ii,jj)

                    IF(ht_sao_base(ii,jj) < 90000. .AND.  &
                       ht_sao_base(ii,jj) < ht_sao_base_min)THEN
                      i_sao = ii
                      j_sao = jj
                      ht_sao_base_min = ht_sao_base(ii,jj)
                    END IF

                  END IF  ! In bounds
                END DO
              END DO

              idelt_index = idelt_index + 1
            END DO
            IF (ht_sao_base_min < -90000.) ht_sao_base_min = ht_sao_base(i,j)

!write(0,'(I2,a,2I4,a,5I4)') myproc,': searching for i/j = ',i,j,' within range ',isrcbgn,isrcend, jsrcbgn, jsrcend, idelt_index
          END IF  ! ht_sao_base > 90000.

          IF (ht_sao_base_min > 90000.) THEN ! Sat. cld but no SAO cld
            n_no_sao2 = n_no_sao2 + 1        ! even in neighbor points.
            cover=sat_cover
            htbase_init = ht_sao_base_min

            IF(t11mu(i,j,1)-t_gnd_k(i,j) < -21.) THEN ! We more likely
                                                      ! have a cloud
              buffer = 2100.
            ELSE
              buffer = sfc_ir_buffer   ! Weed out IR tops < ~5000m AGL
            END IF

!
!-----------------------------------------------------------------------
!
!  Calculate new cloud top and cover
!
!-----------------------------------------------------------------------
!
            cldtop_m_old = cldtop_m(i,j)

            CALL cloud_top(nx,ny,nz,                                    &
                    i,j,zs_1d,t_1d,p_1d,z_ref_lcl,z_lcl(i,j),           &
                    hterain(i,j),tmin,r_missing,                        &
                    t11mu(i,j,2),t_gnd_k(i,j),thresh2,                  &
                    cldtop_m_ir(i,j),cldtop_m(i,j),                     &
                    l_t11mu,l_cloud_present,sat_cover)

!
!-----------------------------------------------------------------------
!
!  Change to cover
!
!-----------------------------------------------------------------------
!
            cover=ir_cover(isat,t11mu(i,j,1),                           &
                              t_gnd_k(i,j),t11mu(i,j,2))
            htbase = MAX(hterain(i,j)+buffer, cldtop_m(i,j)-thk_def )
            i_reset_base = i_reset_base +1

            IF (htbase > cldtop_m(i,j)) THEN
              n_no_sao3 = n_no_sao3 + 1 ! Sat.cld, but no SAO cld,
                                        ! & sat.cld is too low.
            ELSE IF(mod(i_reset_base,50) == 0 .AND. verbose) THEN
              WRITE(6,'(2x,2(a,I4),a)')                                 &
                'Correction at point(',i,', ',j,                        &
                ') because there is low sat.cld, but no SAO cld.'
              WRITE(6,'(2x,2(a,f8.1),2(a,f9.0),a,f5.2)')                &
                'tb=',t11mu(i,j,1),'  tb_cold=',t11mu(i,j,2),           &
                ' ctop=',cldtop_m_old,' ctop_cold=',cldtop_m(i,j),      &
                ' cvr=',cover
            END IF

          ELSE IF (ht_sao_base(i,j) > cldtop_m(i,j)) THEN ! Sat.cld, nut no
                                     ! SAO cld at same pt. Do have
                                     ! SAO cld in nearby pt. AND
                                     ! Sat.top is below ceiling
            cover=sat_cover
            htbase_init = ht_sao_base(i,j)   ! using SAO cld base
            htbase = htbase_init
            cldtop_m_old = cldtop_m(i,j)
            cldtop_m(i,j) = htbase_init + thk_def ! correct sat. cldtop
!
!-----------------------------------------------------------------------
!
!  Find a thinner value for cloud cover consistent with the new
!  higher cloud top and the known brightness temperature.
!  Note that cover is not really used here as an input
!
!-----------------------------------------------------------------------
!
            CALL correct_cover(isat,cover,cover_new,                    &
                     cldtop_m_old,cldtop_m(i,j),i,j,nx,ny,nz,           &
                     zs_1d,t_1d,t11mu(i,j,1),t_gnd_k(i,j),              &
                     istatus)

            IF (istatus /= 1) THEN
              WRITE(6,*)' Error in correct_cover'
              RETURN
            END IF
            cover = cover_new

          ELSE ! Normal use of satellite data
               ! There is a ht_sao_base below cldtop_m in the area
               ! near the grid pt.

            cover=sat_cover
!
!-----------------------------------------------------------------------
!
!  Locate SAO cloud base below satellite cloud top, modify
!  satellite cloud base. Highest SAO ceiling within default
!  thickness range of satellite layer is used.
!
!-----------------------------------------------------------------------
!
            DO k=nz-1,1,-1

              IF (zs(i,j,k) >= htbase_init .AND.                        &
                    zs(i,j,k) <= cldtop_m(i,j)) THEN

                IF (cldcv_sao_wrk(i_sao,j_sao,k)  <= thr_sao_cvr .AND.  &
                    cldcv_sao_wrk(i_sao,j_sao,k+1) > thr_sao_cvr) THEN

!c            We have an SAO base

                  htbase = zs(i,j,k+1)

!c              If SAO (hence satellite) base is above the satellite
!c              cloud top, lower the satellite base by one grid level

                  IF (htbase > cldtop_m(i,j)) htbase = zs(i,j,k)

                  GO TO 301

                END IF ! We have an SAO base
              END IF ! in satellite layer
            END DO ! k
          END IF ! HT_SAO_BASE > 90000.
                 ! No SAO cloud base in the area near the grid pt.

          301  CONTINUE

          IF ( verbose .AND.          &
               htbase /= htbase_init  .AND. mod(i_reset_base,50) == 0 ) THEN
            WRITE(6,'(1x,I4,a)') myproc,': Satellite ceiling reset by SAO.'
            WRITE(6,'(1x,I4,a,2i4,a,f9.0,1x,f9.0,a,f9.0)')              &
                    myproc,': i,j= ',i,j,' htbase new,old=',            &
                    htbase,htbase_init,' cldtop=',cldtop_m(i,j)
          END IF
!
!-----------------------------------------------------------------------
!
!  Add satellite cloud to array
!
!-----------------------------------------------------------------------
!
          DO k=nz,1,-1
            IF (zs(i,j,k) >= htbase .AND.                               &
                zs(i,j,k) <= cldtop_m(i,j)) & ! in satellite layer
            cldcv(i,j,k)=cover
          END DO

        END IF ! l_cloud_present (Cloudy)

      END IF ! non-missing IR data

    END DO ! i=1,nx
  END DO ! j=1,ny
!
  istatus = 1

!
!-----------------------------------------------------------------------
!
!  Write stats on IR insertion
!
!-----------------------------------------------------------------------
!
  IF (verbose) THEN
    CALL mpsumi(n_no_sao1,1)
    CALL mpsumi(n_no_sao2,1)
    CALL mpsumi(n_no_sao3,1)
    IF (myproc == 0) WRITE(6,*)' n_no_sao (1/2/3) = ',n_no_sao1,n_no_sao2,n_no_sao3
  END IF

  CALL compare_rad (nx,ny,nz,clddiag,r_missing,lvldbg,zs,temp_3d,       &
                    isatir,cldcv,t_sfc_k,t_gnd_k,t11mu,                 &
                    cvr_snow,nlyr)
!
!-----------------------------------------------------------------------

  IF (mp_opt > 0) THEN
    CALL mpsendrecv2dew(cldcv,nx,ny,nz,ebc,wbc,1,tem1)
    CALL mpsendrecv2dns(cldcv,nx,ny,nz,nbc,sbc,2,tem1)
  END IF

  DEALLOCATE(zs_wrk, cldcv_sao_wrk)

!if (mp_opt > 0) then
!write(200+myproc,'(6F12.5)') cldcv
!else
!write(200+0,'(6F12.5)') (((cldcv(i,j,k),i=1,35),j=1,35),k=1,nz)
!write(200+1,'(6F12.5)') (((cldcv(i,j,k),i=33,67),j=1,35),k=1,nz)
!write(200+3,'(6F12.5)') (((cldcv(i,j,k),i=1,35),j=33,67),k=1,nz)
!write(200+3,'(6F12.5)') (((cldcv(i,j,k),i=33,67),j=33,67),k=1,nz)
!
!end if

  RETURN
END SUBROUTINE insert_ir

!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE CLOUD_TOP                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cloud_top(nx,ny,nz,                                          &
           i,j,zs_1d,t_1d,p_1d,z_ref_lcl,z_lcl,                         &
           hterain,tmin,r_missing,                                      &
           t11mu,t_gnd_k,thresh2,                                       &
           cldtop_m_ir,cldtop_m,                                       &
           l_t11mu,l_cloud_present,sat_cover)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This routine computes the cloud top height given a 10.7-micron IR
!  brightness temperature and 3D fields of temp and height.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  04/96 Based on the LAPS cloud analysis code (09/95 Dan Birkenheuer)
!
!  MODIFICATION HISTORY:
!
!  04/26/96 (J. Zhang)
!  Modified for the ARPSDAS format.
!  Added full documentation.
!
!  11/01/96 (J. Zhang)
!  Added a scheme for calculation of cloud top height of low clouds
!  (e.g., persistent stratocumulus), where temperature inversion
!  may exist and the simple matching scheme may fail.
!  Reference:
!  Macpherson et al., 1996: The impact of MOPS moisture data in
!      the U.K. Meteorological Office mesoscale data assimilation
!      scheme.   MWR, vol. 124,  1746-1766
!
!  09/10/97 (J. Zhang)
!           Lifting condensation level is input through calling
!           argument.
!
!-----------------------------------------------------------------------
!
!  Variable Declarition
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  For 10.7 micron method
!
!-----------------------------------------------------------------------
!
!  INPUT:
  INTEGER :: nx,ny,nz     ! The model grid size
!
!  INPUT: Vertical 1D arrays
  REAL :: zs_1d(nz)    ! physical height (m) at scalar points
  REAL :: t_1d(nz)  ! 3D temperature analysis on model grid
  REAL :: p_1d(nz)  ! pres. analysis on model grid
!
!  INPUT: At grid point (i,j)
  INTEGER :: i,j          ! grid point location indices
  REAL :: z_lcl           ! Lifting condensation level at i,j
  REAL :: t11mu           ! 10.7 mu IR brightness temperature(K)
  REAL :: t_gnd_k         ! gound temperature
  REAL :: hterain         ! terrain height at the grid point
!
!  INPUT: Constants
  REAL :: tmin
  REAL :: r_missing       ! missing data flag value
  REAL :: thresh2         ! Input from parent routine (=8K)
  REAL :: z_ref_lcl       ! ref. level for computing LCL
!
!  OUTPUT:
  LOGICAL :: l_cloud_present ! "cloudy" indicated
  LOGICAL :: l_t11mu           ! "cloudy" indicated by 10.7-micron IR?
  REAL :: cldtop_m_ir       ! cld top height(m) obtained from t11mu
  REAL :: cldtop_m           ! cld top height(m) obtained from IR data
  REAL :: sat_cover          ! cloud fraction obtained from IR data
!
!  LOCAL:
  REAL :: gamma_s         ! ref. moist adiabatic lapse rate
!
!  LOCAL:  For the conceptual model for the cloud top
  REAL :: t_base,p_base_pa,zbase,ztop,dzinc,zht
  PARAMETER(dzinc=100.0)
  REAL :: p,press,temp,eso,dlvdt,des,dtz,rl,es
  REAL :: t11mu_loc
  INTEGER :: nlevel,nlm1
!
!  CONSTANTS:
  REAL :: gamma_d   ! dry adiabatic lapse rate (K/m)
  PARAMETER (gamma_d = 9.8/1004.0)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: kl,j1
  REAL :: arg,frac_k,frac_z,z_ref,t_ref
!
!-----------------------------------------------------------------------
!
!  Include files.
!
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Define constants
!
!-----------------------------------------------------------------------
!
  dlvdt=-2.3693E+3        ! J/kg/K
  eso=610.78              ! pa
!
!-----------------------------------------------------------------------
!
!  This section finds the cloud top using 10 micorn IR data
!
!-----------------------------------------------------------------------
!
  cldtop_m_ir = r_missing ! or zeros
  l_t11mu = .false.
!
  IF (t11mu - t_gnd_k < -thresh2) THEN ! probably clds
    l_t11mu = .true.
    IF (t11mu < 253.15) GO TO 300 ! Using scan&match method
!
!-----------------------------------------------------------------------
!
!  When Tb >= -20 deg C, find the cloud top height using a
!  conceptual model presented in Macpherson et al. (1996).
!
!-----------------------------------------------------------------------
!
    zbase = z_lcl
    IF (zbase >= zs_1d(nz-1)) GO TO 300
! Using scanning & matching method
!
!-----------------------------------------------------------------------
!
!  Find the model pressure and temperature at the lifting
!  condensation level.
!
!-----------------------------------------------------------------------
!
    t_ref=-999.
    DO kl = 2,nz-1
      z_ref = z_ref_lcl + hterain
      IF (z_ref <= zs_1d(kl)) THEN
        frac_z = (z_ref-zs_1d(kl-1))/(zs_1d(kl)-zs_1d(kl-1))
        t_ref = t_1d(kl-1) + frac_z * (t_1d(kl) - t_1d(kl-1))
        EXIT
      END IF
    END DO  ! kl
    IF(t_ref < -900.) THEN
      WRITE(6,'(a)') ' Error for Tbase,aborting....'
      STOP
    END IF
    t_base = t_ref - (z_lcl -z_ref)*gamma_d
!
    p_base_pa=-999.
    DO kl = 2,nz-1
      IF (zbase <= zs_1d(kl)) THEN
        frac_z = (zbase-zs_1d(kl-1))/(zs_1d(kl)-zs_1d(kl-1))
        arg = ALOG(p_1d(kl-1))                                          &
               + frac_z * (ALOG(p_1d(kl))-ALOG(p_1d(kl-1)))
        p_base_pa = EXP(arg)
        EXIT
      END IF
    END DO  ! kl
    IF( p_base_pa < 0. ) THEN
      PRINT*,'Error for Zbase,aborting....'
      STOP
    END IF

    IF (t_base <= t11mu) GO TO 300
!
!-----------------------------------------------------------------------
!
!  Using scanning & matching method
!
!  Find the cloud top height, which is the level where the temp.
!  reaches the satellite brightness temperature throught the moist
!  adiabatic cooling process.
!
!-----------------------------------------------------------------------
!
    temp   = t_base    !  temp   = t_base_k
    press  = p_base_pa
    gamma_s = 0.004        ! K/m
    nlevel = (temp -t11mu)/gamma_s/dzinc  ! a guess for cloud top
    IF(nlevel < 0) THEN
      WRITE(6,'(a)') ' Error in zlevel. aborting....'
      WRITE(6,'(a,i6,a,f9.1,a,f9.3,/a,f9.1,a,f9.2)')                    &
        ' nlevel=',nlevel,' temp=',temp,' gam=',gamma_s,                &
             ' dzinc=',dzinc,' IR temp=',t11mu
      STOP
    END IF
    nlm1   = nlevel
    IF(nlm1 < 1) nlm1=1
    zht    = zbase
!
    DO j1=1,nlm1
      rl   = lathv+(273.15-temp)*dlvdt    ! Lv as func of Temp.
      arg  = rl*(temp-273.15)/273.15/temp/rv
      es   = eso*EXP(arg)                 ! satur. vapor press.
      arg  = -g*dzinc/rd/temp
      p    = press*EXP(arg)               ! hydrosta. press.decrease
!
!-----------------------------------------------------------------------
!
!  Calculate saturated adiabatic lapse rate
!
!-----------------------------------------------------------------------
!
      des   = es*rl/temp/temp/rv
      dtz   = -g*((1.0+0.621*es*rl/(press*rd*temp))/                    &
               (cp+0.621*rl*des/press))
      IF(ABS(dtz) < 1.0E-15) THEN
        PRINT*,' Zero dt/dz:',dtz,' g=',g,' es=',es,' press=',press
        PRINT*,'temp=',temp,' cp=',cp,' rl=',rl,' des=',des
        STOP
      END IF
      zht   = zht+dzinc  ! moist adiabatic ascent every 100m
      temp  = temp+dtz*dzinc
      IF(temp <= t11mu) THEN   ! Cloud top is found
        ztop = zht + (t11mu-temp)/dtz
        GO TO 150
      END IF
      press = p
    END DO  ! j=1,nlm1
    ztop = zht + (t11mu-temp)/dtz
    150     CONTINUE
!
!-----------------------------------------------------------------------
!
!  Apply quality check to the cloud top height.
!
!-----------------------------------------------------------------------
!
    IF(ztop <= zs_1d(1)) THEN
      PRINT*,' Error, ztop for cloud top model is out of'               &
             ,' bound at grid pt(',i,j,').'
      PRINT*,' ztop=',ztop,' hterain=',hterain,                               &
             ' z_bottom=',zs_1d(1)
      PRINT*,'aborting...'
      STOP
    END IF
!
    IF(ztop > zs_1d(nz-1)) THEN
      PRINT*,' Warning, ztop for cloud top model is out of'             &
               ,' bound at grid pt(',i,j,').'
      PRINT*,' ztop=',ztop,' hterain=',hterain                                &
               ,' z_top=',zs_1d(nz-1)
      PRINT*,' ztop is set to the highest model level'
      cldtop_m_ir = zs_1d(nz-1)
      GO TO 900
    END IF
    cldtop_m_ir = ztop
    GO TO 900
!
    300     CONTINUE
!
!-----------------------------------------------------------------------
!
!  When Tb < -20 deg C, use the scheme that finds the cloud top
!  height by scanning the model temperature profile for a value
!  matching IR temperature.
!
!-----------------------------------------------------------------------
!
    t11mu_loc=t11mu
    IF(t11mu < tmin) THEN
      WRITE(6,'(a)') ' Warning: cloud top colder than tmin in column'
      WRITE(6,'(a)') '   i  j   t11mu   t_top'
      WRITE(6,'(1x,2i3,2f8.2)') i,j,t11mu,t_1d(nz-1)
      WRITE(6,'(a)') ' Cloud top temperature is set to tmin'
      t11mu_loc=tmin
    END IF
!
    DO kl = nz-2,1,-1
!
!-----------------------------------------------------------------------
!
!  Find the lowest temp. crossing point in the model temperature
!  profile.
!
!-----------------------------------------------------------------------
!
      IF( (t_1d(kl)-t11mu_loc) *  (t_1d(kl+1)-t11mu_loc) <= 0.) THEN   ! Crossing Pt

        frac_k = (t11mu_loc - t_1d(kl))                                   &
                    / (t_1d(kl+1) - t_1d(kl))
        arg = zs_1d(kl) + frac_k * (zs_1d(kl+1)-zs_1d(kl))

        IF(arg >= hterain) THEN
          cldtop_m_ir = arg
        ELSE
          WRITE(6,*)' Warning: Cloud Top Below Ground - flagged'

          WRITE(6,'(a,a)') '  i  j kl  t11mu_loc  t_abv  t_blw frac',     &
           '  hgt_m    hterain  ctop_8 '
          WRITE(6,'(1x,3i3,3f7.2,f5.2,3f8.0)') i,j,kl,t11mu_loc,          &
            t_1d(kl+1),t_1d(kl),frac_k,arg,hterain,cldtop_m_ir
        END IF

      END IF   ! Crossing point found
    END DO ! kl

  ELSE ! No clouds according to IR satellite data

    l_t11mu = .false.

  END IF   ! t11mu-t_gnd_k .lt. -thresh2
!
!-----------------------------------------------------------------------
!
!  Set variables
!
!-----------------------------------------------------------------------
!
  900   CONTINUE
  l_cloud_present = l_t11mu
  cldtop_m = cldtop_m_ir
  sat_cover = 1.0

  RETURN
END SUBROUTINE cloud_top

!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE CORRECT_COVER                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE correct_cover(isat,cover_in,cover_new_f,                    &
                         cldtop_old,cldtop_new,i,j,nx,ny,nz,           &
                         zs_1d,t_1d,t11mu,t_gnd_k,                     &
                         istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Find a thinner value for cloud cover consistent with the new
!  higher cloud top and the known brightness temperature.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Dan Birkenheuer?)
!  09/95.
!
!  MODIFICATION HISTORY:
!
!  04/29/96 (J. Zhang)
!  Modified for the ARPSDAS format.
!  Added full documentation.
!
!-----------------------------------------------------------------------
!
!  Variable Declaration
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!  INPUT:
  INTEGER :: isat
  INTEGER :: i,j         ! grid point index
  INTEGER :: nx,ny,nz    ! grid size
  REAL :: t_1d(nz)       ! temp. analysis field
  REAL :: zs_1d(nz)      ! phycisal heights (m) for a column of scalar
                         ! grid points

  REAL :: t11mu          ! satellite 10.7 micron IR brightness temp (K)
  REAL :: t_gnd_k        ! ground skin temp.
  REAL :: cover_in       ! cloud cover before correction
  REAL :: cldtop_old     ! cloud top height before correction
!
!  OUTPUT:
  INTEGER :: istatus
  REAL :: cover_new_f    ! cloud cover after correction
  REAL :: cldtop_new     ! cloud top height after correction
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL :: z_temp,temp_old,temp_new,frac
  REAL :: ir_cover    ! a function
  REAL :: cover_old,cover_new
  INTEGER :: iz_temp
  INTEGER :: iwrite
  DATA iwrite /0/
  SAVE iwrite

  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Find Temperature of old cloud top
!
!-----------------------------------------------------------------------
!
  CALL hgt_to_zk (cldtop_old,z_temp,nz,zs_1d,istatus)
  IF(istatus /= 1 .AND. myproc == 0)THEN
    PRINT*,' Old cldtop is out of domain range:', z_temp
    PRINT*,' Old cldtop=',cldtop_old,' zs range from',zs_1d(1)          &
            ,' to ',zs_1d(nz-1)
  END IF

  z_temp = MAX(1.,MIN(z_temp,FLOAT(nz-1)-.001))
  iz_temp = INT(z_temp)
  frac = z_temp - iz_temp
  temp_old = t_1d(iz_temp)    * (1. - frac)                             &
          +  t_1d(iz_temp+1)  * frac
!
!-----------------------------------------------------------------------
!
!  Find Temperature of new cloud top
!
!-----------------------------------------------------------------------
!
  CALL hgt_to_zk (cldtop_new,z_temp,nz,zs_1d,istatus)
  IF(istatus /= 1 .AND. myproc == 0)THEN
    PRINT*,' New cldtop is out of domain range:', z_temp
    PRINT*,' New cldtop=',cldtop_new,' zs range from',zs_1d(1)          &
            ,' to ',zs_1d(nz-1)
!        return
  END IF

  z_temp = MAX(1.,MIN(z_temp,FLOAT(nz-1)-.001))
  iz_temp = INT(z_temp)
  frac = z_temp - iz_temp
  temp_new = t_1d(iz_temp) * (1. - frac)                                &
              + t_1d(iz_temp+1) * frac
!
!-----------------------------------------------------------------------
!
!  This one utilizes a linear approximation to
!  the sigma T**4 relationship
!
!-----------------------------------------------------------------------
!
  cover_old = MIN(cover_in,1.0)
  cover_new = cover_old * (t11mu-t_gnd_k)/(temp_new-t_gnd_k)
!
!-----------------------------------------------------------------------
!
!  This one utilizes the sigma T**4 relationship
!
!-----------------------------------------------------------------------
!
  cover_new_f = ir_cover(isat,t11mu,t_gnd_k,temp_new)

  IF((j-1) == INT(j-1)/10*10) THEN
    iwrite = iwrite + 1
    IF(iwrite < 15 .AND. myproc == 0) THEN
      WRITE(6,'(1x,a,a)') '**CORR_CVR**  i  j  tg_k   told   tnew ',    &
                          '  ctold   ctnew   cvnew  cvrnewf'
      WRITE(6,'(1X,12X,2I3,3F7.0,2F8.0,2F8.2)') i,j,t_gnd_k,            &
        temp_old,temp_new,cldtop_old,cldtop_new,cover_new,cover_new_f
    END IF
  END IF

  istatus=1

  RETURN
END SUBROUTINE correct_cover

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION IR_COVER                    ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION ir_cover(isat,t11mu,t_gnd_k,t_cld)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Estimate cloud cover based on IR temperature compared to
!  ground temperature.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  07/95
!
!  MODIFICATION HISTORY:
!  05/01/96  (Jian Zhang)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!  INPUT:
!
  INTEGER, INTENT(IN) :: isat
  REAL, INTENT(IN) :: t11mu       ! satellite 10.7 micron bright temp. (K)
  REAL, INTENT(IN) :: t_gnd_k     ! ground skin temp. (K)
  REAL, INTENT(IN) :: t_cld       ! estimated brightness temp. with SAO cld
!
!  OUTPUT:
!
  REAL :: ir_cover ! modified cld cover using !0-micron data
!
!  LOCAL:
!
  REAL :: r_sfc,r_sat,r_cld
  REAL :: temp_to_rad       ! a function
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  r_sfc = temp_to_rad(isat,t_gnd_k)
  r_sat = temp_to_rad(isat,t11mu)
  r_cld = temp_to_rad(isat,t_cld)

  ir_cover = (r_sat-r_sfc) / (r_cld-r_sfc)

  RETURN
  END FUNCTION ir_cover

!
!##################################################################
!##################################################################
!######                                                      ######
!######              FUNCTION TEMP_TO_RAD                    ######
!######                                                      ######
!##################################################################
!##################################################################
!
  FUNCTION temp_to_rad(isat,temp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Convert the radiance to the brightness temperature
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  07/95
!
!  MODIFICATION HISTORY:
!  05/01/96  (Jian Zhang)
!  03/11/03  (Keith Brewster)
!            Re-written to use latest version of calibration
!            constants from CIMSS and streamlined.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: temp_to_rad
  INTEGER, INTENT(IN) :: isat
  REAL, INTENT(IN)  :: temp
!
!-----------------------------------------------------------------------
!
! Include file for calibration common block.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'adassat.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  temp_to_rad=fk1(isat)/                                                &
              (exp(fk2(isat)/(bc1(isat)+(bc2(isat)*temp)))-1.0)

  RETURN
  END FUNCTION temp_to_rad
!
!##################################################################
!##################################################################
!######                                                      ######
!######              FUNCTION RAD_TO_TEMP                    ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION rad_to_temp(isat,rad)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Convert the radiance to the brightness temperature
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  07/95
!
!  MODIFICATION HISTORY:
!  05/01/96  (Jian Zhang)
!  03/11/03  (Keith Brewster)
!            Re-written to use latest version of calibration
!            constants from CIMSS, to use multiple satellites,
!            and streamlined.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: rad_to_temp
  INTEGER, INTENT(IN) :: isat
  REAL, INTENT(IN)  :: rad
!
!-----------------------------------------------------------------------
!
! Include file for calibration common block.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'adassat.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  rad_to_temp=(fk2(isat)/(alog((fk1(isat)/rad)+1.0))-bc1(isat))         &
              *bc2inv(isat)
!
  RETURN
END FUNCTION rad_to_temp

!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE COMPARE_RAD                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE compare_rad (nx,ny,nz,clddiag,r_missing,lvldbg,zs,t_3d,      &
                        isatir,cldcv,t_sfc_k,t_gnd_k,t11mu,             &
                        cvr_snow,nlyr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This routine compares the analyzed clouds to the 10.7 micron radiation
!  and determines adjusted cloud fractions of the cloud layers to yield
!  a better fit.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (?)
!  09/95.
!
!  MODIFICATION HISTORY:
!
!  05/01/96 J. Zhang
!           Modified for the ARPSDAS format.
!           Added full documentation.
!  05/01/99 K. Brewster
!           Added check for missing IR data.
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
  INTEGER :: nx,ny,nz     ! grid size
!
!  INPUT:
  INTEGER :: clddiag      ! print diagnostic fields?
  REAL :: r_missing       ! bad flag
  REAL :: t_3d(nx,ny,nz)  ! temp. analysis on the model grid
  REAL :: zs(nx,ny,nz)    ! physical heights (m) at model scalar pts.

  INTEGER :: isatir(nx,ny)
  REAL :: cvr_snow(nx,ny) ! sfc snow coverage fraction
  REAL :: t11mu(nx,ny)    ! satellite band8 brightness temp.
  REAL :: t_gnd_k(nx,ny)  ! ground skin temp.
  REAL :: t_sfc_k(nx,ny)  ! surface air temperature
!
!  OUTPUT:
  REAL :: cldcv(nx,ny,nz)   ! adjusted cloud cover analysis
!
!  LOCAL:
  REAL :: zs_1d(nz)             ! phy. heights (m) in a column of grid
  REAL :: t_1d(nz)              ! temp in a column of grid
  REAL :: cldcv_1d(nz)      ! cloud cover in one culomn of grid.
  REAL :: cldcv_1dref(nz)   ! cloud cover in one culomn of grid.

  REAL :: a(nz)                ! max. cldcvr in each cloud layer
  REAL :: f(nz)
  REAL :: a_new(nz)
  REAL :: f_new(nz)

  INTEGER :: ilyr(nz)       ! Dims needs to be changed to nz
  INTEGER :: ilyr_new(nz)   ! Dims needs to be changed to nz

  INTEGER :: nlyr,nlyr_new
  LOGICAL :: l_correct

  INTEGER, INTENT(IN) :: lvldbg
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL :: corr_thr,cover_step
  REAL :: delta_cover,delta_cover_ref
  INTEGER :: i_correct,iter,iter_max,n_iter
  REAL :: t_effective,tdiff,tdiff_ref
  INTEGER :: i,j,k,l,iwrite,imid,jmid,kmid,ibeg,iend,isat
  INTEGER :: n_clear,n_clear_ns,n_clear_sn,n_cloudy,n_total             &
          ,n_correct
  REAL :: tdiff_sumsq,tdiff_cld_sumsq,tdiff_corr_sumsq                  &
      ,tdiff_corr_sumsq2,tdiff_corr_sumsq3,tdiff_corr_cld_sumsq
  REAL :: tir_g_clr_sum,tir_a_clr_sum,tir_g_clr_sumsq                   &
      ,tir_a_clr_sumsq,tir_g_clr_ns_sum,tir_a_clr_ns_sum                &
      ,tir_g_clr_ns_sumsq,tir_a_clr_ns_sumsq,tir_g_clr_sn_sum           &
      ,tir_a_clr_sn_sum,tir_g_clr_sn_sumsq,tir_a_clr_sn_sumsq

  REAL :: frac,frac_clouds

  LOGICAL :: verbose = .FALSE.

  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (lvldbg > 90 .AND. myproc == 0) verbose = .TRUE.

  imid=nx/2
  jmid=ny/2
  ibeg=MAX( 1,imid-5)
  iend=MIN(nx,imid+5)
  kmid=nz/2
  IF(clddiag == 1) THEN
    WRITE(6,'(a)')' Comparing radiation'
    WRITE(6,611) nx,ny,nz
    611     FORMAT(1X,' nx=',i3,' ny=',i3,' nz=',i3)
    WRITE(6,'(a)')'grid point heights:'
    WRITE(6,613) (zs(i,jmid,kmid),i=ibeg,iend)
    613     FORMAT(10F8.0)
    WRITE(6,'(a)')'grid point temp:'
    WRITE(6,613) (t_3d(i,jmid,kmid),i=ibeg,iend)
    WRITE(6,'(a)')
!
    WRITE(6,'(a)') 'Surface temperatures: j=jmid, i=imid-5,+5'
    WRITE(6,614) (t_sfc_k(i,jmid),i=ibeg,iend)
    614     FORMAT(1X,10F7.0)
    WRITE(6,'(a)') 'Ground skin temperatures: j=jmid, i=imid-5,+5'
    WRITE(6,614) (t_gnd_k(i,jmid),i=ibeg,iend)
    WRITE(6,'(a)') 'IR sat temperatures: j=jmid, i=imid-5,+5'
    WRITE(6,614) (t11mu(i,jmid),i=ibeg,iend)
    WRITE(6,'(a)')
  END IF

  corr_thr = 4.             ! if diff of t11mu-t11mue > 4K, correction
  cover_step = 0.05         ! add/minus 0.05 cldcvr each correction
  iter_max = 10             ! we'll do maximum of 10 iterns

  iwrite = 0                ! output message flag
  n_clear = 0               !
  n_clear_ns = 0
  n_clear_sn = 0
  n_cloudy = 0
  n_total = 0
  n_correct = 0
  n_iter = 0
  tdiff_sumsq = 0.          ! sum of sqr(tIRo-t11mue) for all pts
                             ! before correction
  tdiff_cld_sumsq = 0.      ! sum of sqare(tIRo-t11mue) for cloudy pts
                             ! before correction
  tdiff_corr_sumsq3 = 0.    ! sum of sqr(tIRo-t11mue) for cloudy pts
                            ! to be corrected, before the correction
  tdiff_corr_sumsq = 0.     ! sum of sqr(tIRo-t11mue) for all pts
                             ! after correction
  tdiff_corr_cld_sumsq = 0. ! sum of sqr(tIRo-t11mue) for all cldy pts
                             ! after correction
  tdiff_corr_sumsq2 = 0.    ! sum of sqr(tIRo-t11mue) for cloudy pts
                            ! to be corrected, after the correction
  tir_g_clr_sum   = 0.      ! sum of tir-tg for clear air pts
  tir_a_clr_sum   = 0.      ! sum of tir-tsfc for clear air pts
  tir_g_clr_sumsq = 0.      ! sum of sqare(tir-tg) for clear air pts
  tir_a_clr_sumsq = 0.      ! sum of sqare(tir-tsfc) for clear air pts
  tir_g_clr_ns_sum   = 0.   ! for clear air with < 0.5 snow cover pts
  tir_a_clr_ns_sum   = 0.
  tir_g_clr_ns_sumsq = 0.
  tir_a_clr_ns_sumsq = 0.
  tir_g_clr_sn_sum   = 0.   ! for clear air with > 0.5 snow cover pts
  tir_a_clr_sn_sum   = 0.
  tir_g_clr_sn_sumsq = 0.
  tir_a_clr_sn_sumsq = 0.

!
!-----------------------------------------------------------------------
!
!  Compare and correct (if necessary) the cloud cover analysis field
!  for each grid point.
!
!-----------------------------------------------------------------------
!
  DO j = 1,ny-1
    DO i = 1,nx-1
      isat=isatir(i,j)
      IF(t11mu(i,j) > 0.) THEN
        DO k = 1,nz
          zs_1d(k) = zs(i,j,k)
          t_1d(k) = t_3d(i,j,k)
        END DO

        DO k = 1,nz
          cldcv_1dref(k) = cldcv(i,j,k)
        END DO

        CALL cvr_to_t11mue (nz,isat,zs_1d,t_1d,                           &
                   cldcv_1dref,t_gnd_k(i,j),                            &
                   a,f,ilyr,t_effective,nlyr)
!
        tdiff = t11mu(i,j)-t_effective ! IR sat temp - calculated
        frac_clouds = 1.-f(nlyr)   ! column total cloud cover

        IF(nlyr == 1) THEN
          tir_g_clr_sum   = tir_g_clr_sum   + tdiff
          tir_g_clr_sumsq = tir_g_clr_sumsq + tdiff**2
          tir_a_clr_sum   = tir_a_clr_sum + t11mu(i,j)-t_sfc_k(i,j)
          tir_a_clr_sumsq=tir_a_clr_sumsq+(t11mu(i,j)-t_sfc_k(i,j))**2
          n_clear = n_clear + 1

          IF(cvr_snow(i,j) /= r_missing) THEN
            IF(cvr_snow(i,j) < 0.5) THEN
              tir_g_clr_ns_sum   = tir_g_clr_ns_sum   + tdiff
              tir_g_clr_ns_sumsq = tir_g_clr_ns_sumsq + tdiff**2
              tir_a_clr_ns_sum = tir_a_clr_ns_sum                       &
                                  + t11mu(i,j)-t_sfc_k(i,j)
              tir_a_clr_ns_sumsq = tir_a_clr_ns_sumsq                   &
                                  + (t11mu(i,j) - t_sfc_k(i,j))**2
              n_clear_ns = n_clear_ns + 1
            ELSE
              tir_g_clr_sn_sum   = tir_g_clr_sn_sum   + tdiff
              tir_g_clr_sn_sumsq = tir_g_clr_sn_sumsq + tdiff**2
              tir_a_clr_sn_sum = tir_a_clr_sn_sum                       &
                                  + t11mu(i,j)-t_sfc_k(i,j)
              tir_a_clr_sn_sumsq = tir_a_clr_sn_sumsq                   &
                                  +(t11mu(i,j)-t_sfc_k(i,j))**2
              n_clear_sn = n_clear_sn + 1
            END IF
          END IF
        ELSE
          n_cloudy = n_cloudy + 1
          tdiff_cld_sumsq = tdiff_cld_sumsq + tdiff**2
        END IF

        n_total = n_total + 1
        tdiff_sumsq = tdiff_sumsq + tdiff**2
!
!-----------------------------------------------------------------------
!
!    Apply corrections?
!
!    Only when the following four conditions are satisfied:
!    1. at least one cloud layer
!    2. the difference between the observed and estimated tb8
!       is larger than the threshold (4K)
!    3. tirobs - t_gnd < -8K
!    4. column total cloud cover is larger than 40%
!
!-----------------------------------------------------------------------
!
        IF (nlyr >= 2 .AND. ABS(tdiff) > corr_thr                       &
              .AND. t_gnd_k(i,j) - t11mu(i,j) > 8.                      &
              .AND. frac_clouds > 0.4 )  THEN
          l_correct = .true.
          n_correct = n_correct + 1
          tdiff_corr_sumsq3 = tdiff_corr_sumsq3 + tdiff**2
        ELSE
          l_correct = .false.
        END IF
!
!-----------------------------------------------------------------------
!
!    This corrective section is now turned on
!
!-----------------------------------------------------------------------
!
        iter = 0
        delta_cover = 0.

        905 IF (l_correct) THEN
          iter = iter + 1
          tdiff_ref = tdiff             ! save initial difference
          delta_cover_ref = delta_cover ! save initial correction

          IF(tdiff < 0.) THEN
            delta_cover = delta_cover + cover_step
          ELSE
            delta_cover = delta_cover - cover_step
          END IF

          CALL apply_correction (cldcv_1dref,cldcv_1d,                  &
                   nz,ilyr,f,delta_cover)

          CALL cvr_to_t11mue (nz,isat,zs_1d,t_1d,                         &
                   cldcv_1d,t_gnd_k(i,j),                               &
                   a_new,f_new,ilyr_new,t_effective,nlyr_new)
!
          tdiff = t11mu(i,j)-t_effective
          frac_clouds = 1.-f_new(nlyr_new)
!
!-----------------------------------------------------------------------
!
!  Continue to apply corrections?
!
!-----------------------------------------------------------------------
!
          IF(nlyr_new >= 2 .AND. frac_clouds > 0.4) THEN
            i_correct = 1
          ELSE
            i_correct = 0
          END IF

          IF(MOD(iwrite,50) == 0 .AND. verbose) THEN
            iwrite = iwrite + 1
            WRITE(6,*) 'iter=',iter ,' tdiff_ref=',tdiff_ref            &
                       ,'tg=',t_gnd_k(i,j)
            WRITE(6,*) (a(l),l=nlyr-1,1,-1)
            WRITE(6,641)
            641 FORMAT(1X,' i  j  ncn  tb8o   t11mue  tdiff  cldcvm i_c'&
                       ,6(' cvrln'))
            WRITE(6,642,ERR=912)i,j,nlyr_new-1,t11mu(i,j),t_effective,  &
                                tdiff,frac_clouds,i_correct,            &
                                (a_new(l),l=nlyr_new-1,1,-1)
            642 FORMAT(1X,2I3,i4,2F7.0,f7.1,f7.3,i2,2X,10F6.2)
            912 CONTINUE
          END IF

          IF(i_correct == 1 .AND. iter < iter_max                       &
              .AND. ABS(tdiff) < ABS(tdiff_ref)                         &
              .AND. tdiff * tdiff_ref > 0.) GO TO 905
! Loop back & increment cover
          n_iter = n_iter + iter
!
!-----------------------------------------------------------------------
!
!      Final iteration
!
!-----------------------------------------------------------------------
!
          IF(.false. .OR. tdiff*tdiff_ref >= 0.) THEN
                                      ! Select best of last two increments
            IF(ABS(tdiff) >= ABS(tdiff_ref)) THEN
              delta_cover = delta_cover_ref
              CALL apply_correction (cldcv_1dref,cldcv_1d               &
                    ,nz,ilyr,f,delta_cover)
              tdiff = tdiff_ref
            END IF
          ELSE           ! Do one Newton iteration
            frac = tdiff / (tdiff - tdiff_ref)
            delta_cover=delta_cover+frac*(delta_cover_ref-delta_cover)

            CALL apply_correction (cldcv_1dref,cldcv_1d                 &
                  ,nz,ilyr,f,delta_cover)

            CALL cvr_to_t11mue (nz,isat,zs_1d,t_1d,                       &
                   cldcv_1d,t_gnd_k(i,j),                               &
                   a_new,f_new,ilyr_new,t_effective,nlyr_new)

            tdiff = t11mu(i,j)-t_effective
            n_iter = n_iter + 1
          END IF
!
!-----------------------------------------------------------------------
!
!  Put the corrected column cloud coverback to 3D cloud cover field
!
!-----------------------------------------------------------------------
!
          DO k=1,nz
            cldcv(i,j,k) = cldcv_1d(k)
          END DO

        END IF ! Corrections were made

        tdiff_corr_sumsq     = tdiff_corr_sumsq + tdiff*tdiff
        IF(nlyr > 1)      & ! cloudy (at least before adjustment)
        tdiff_corr_cld_sumsq = tdiff_corr_cld_sumsq + tdiff**2

        IF(l_correct) THEN
          tdiff_corr_sumsq2 = tdiff_corr_sumsq2 + tdiff*tdiff
        END IF

      END IF

    END DO ! I
  END DO ! J
!
!-----------------------------------------------------------------------
!
!  Write out statistics on consistency of IR sat data and cloud cover
!
!-----------------------------------------------------------------------
!
  IF(n_clear > 0 .AND. verbose) THEN
    WRITE(6,'(/a,i5,2f9.3)')                                            &
     ' Mean/RMS IR sat/gnd temp residual in clear skies =',             &
       n_clear,tir_g_clr_sum/FLOAT(n_clear),                            &
       SQRT(tir_g_clr_sumsq/FLOAT(n_clear))
    WRITE(6,'(a,i5,2f9.3)')                                             &
     ' Mean/RMS IR sat/air temp residual in clear skies =',             &
       n_clear,tir_a_clr_sum/FLOAT(n_clear),                            &
       SQRT(tir_a_clr_sumsq/FLOAT(n_clear))
  END IF

  IF(n_clear_ns > 0 .AND. verbose) THEN
    WRITE(6,'(/a,i5,2f9.3)')                                            &
     ' Mean/RMS IR sat/gnd temp resid in clear/nsnow skies =',          &
       n_clear_ns,tir_g_clr_ns_sum/FLOAT(n_clear_ns),                   &
       SQRT(tir_g_clr_ns_sumsq/FLOAT(n_clear_ns))
    WRITE(6,'(a,i5,2f9.3)')                                             &
     ' Mean/RMS IR sat/air temp resid in clear/nsnow skies =',          &
       n_clear_ns,tir_a_clr_ns_sum/FLOAT(n_clear_ns),                   &
       SQRT(tir_a_clr_ns_sumsq/FLOAT(n_clear_ns))
  END IF

  IF(n_clear_sn > 0 .AND. verbose) THEN
    WRITE(6,'(/a,i5,2f9.3)')                                            &
     ' Mean/RMS IR sat/gnd temp resid in clear/snow skies =',           &
       n_clear_sn,tir_g_clr_sn_sum/FLOAT(n_clear_sn)
    WRITE(6,'(a,i5,2f9.3)')                                             &
     ' Mean/RMS IR sat/air temp resid in clear/snow skies =',           &
       n_clear_sn,tir_a_clr_sn_sum/FLOAT(n_clear_sn),                   &
       SQRT(tir_a_clr_sn_sumsq/FLOAT(n_clear_sn))
  END IF

  IF(n_total > 0 .AND. verbose) THEN
    WRITE(6,'(/a,i5,2f9.3)')                                            &
     ' RMS IR sat/teff residual (bfr corr) in all skies =',             &
     n_total,SQRT(tdiff_sumsq/FLOAT(n_total))
    WRITE(6,'(a,i5,f9.3)')                                              &
     ' RMS IR sat/teff residual (aft corr) in all skies =',             &
     n_total,SQRT(tdiff_corr_sumsq/FLOAT(n_total))
  END IF

  IF(n_cloudy > 0 .AND. verbose) THEN
    WRITE(6,'(/a,i5,f9.3)')                                             &
     ' RMS IR sat/teff residual (bfr corr) in cldy skies =',            &
     n_cloudy,SQRT(tdiff_cld_sumsq/FLOAT(n_cloudy))
    WRITE(6,'(a,i5,f9.3)')                                              &
     ' RMS IR sat/teff residual (aft corr) in cldy skies =',            &
     n_cloudy,SQRT(tdiff_corr_cld_sumsq/FLOAT(n_cloudy))
  END IF

  IF(n_correct > 0 .AND. verbose) THEN
    WRITE(6,'(/a,i5,f9.3)')                                             &
     ' RMS IR sat/teff resid (bfr corr - corrected pts) =',             &
     n_correct,SQRT(tdiff_corr_sumsq3/FLOAT(n_correct))
    WRITE(6,'(a,i5,f9.3)')                                              &
     ' RMS IR sat/teff resid (aft corr - corrected pts) =',             &
     n_correct,SQRT(tdiff_corr_sumsq2/FLOAT(n_correct))
    WRITE(6,'(a,i12,f6.2)') ' Total/Average # of iterations = ',        &
     n_iter,FLOAT(n_iter)/FLOAT(n_correct)
  END IF

  RETURN
END SUBROUTINE compare_rad

!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE CVR_TO_T10MUE               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cvr_to_t11mue (nz,isat,zs_1d,t_1d,                           &
                        cvr,t_gnd_k,                                    &
                        a,f,ilyr,t_effective,nlyr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  To calculate the estimated band8 brightness temperature using
!  the cloud cover analysis.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  (Jian Zhang)
!  04/1996   Based on the LAPS cloud analysis code of 09/95
!
!  MODIFICATION HISTORY:
!
!  04/26/96  J. Zhang
!            Modified for the ARPSDAS format. Added full
!            documentation.
!  07/27/96  J. Zhang
!            Modified code so that the cloud layer with top at the
!            upper boundary of the cloud grid is counted.
!  04/11/97  J. Zhang
!            Include adascld24.inc for ncloud
!  08/06/97  J. Zhang
!            Change adascld24.inc to adascld25.inc
!  09/11/97  J. Zhang
!            Change adascld25.inc to adascld26.inc
!  04/20/98  J. Zhang
!            Based on the arps4.3.3 version. Abandoned cloud grid,
!            using the model grid.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!  INPUT
!
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: isat
  REAL, INTENT(IN) :: cvr(nz)     ! Cloud cover from analysis

  REAL, INTENT(IN) :: zs_1d(nz)
  REAL, INTENT(IN) :: t_1d(nz)
  REAL, INTENT(IN) :: t_gnd_k
!
!  OUTPUT:
!
  REAL, INTENT(OUT) :: t_effective  ! effective brightness temp. estimated
                                    ! from cloud cover analysis
  INTEGER, INTENT(OUT) :: nlyr
  INTEGER, INTENT(OUT) :: ilyr(nz)  ! Layer index for each cloud lvl (needs nz)

  REAL, INTENT(OUT) :: a(nz)   ! Cloud fractions of layers
  REAL, INTENT(OUT) :: f(nz)   ! Apparnt "x-sectn" of cldlyrs seen from above
!
!  LOCAL:
!
  INTEGER :: ik(nz)      ! Height level representative of cloud layers
  REAL :: temp_lyr(nz)   ! temp. at "ik" lavels
!
!  FUNCTIONS:
!
  REAL :: temp_to_rad
  REAL :: rad_to_temp
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k,n
  REAL :: rsum,sumf
  REAL :: cvr_thresh
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
  nlyr = 0

  IF(cvr(nz) >= cvr_thresh) THEN  ! cldtop at the upper boundary
    nlyr = nlyr + 1
    a(nlyr) = cvr(nz)
    ik(nlyr) = nz
    ilyr(nz) = nlyr
  ELSE
    ilyr(nz)=0
  END IF

  DO k = nz-1,1,-1
    IF(cvr(k) >= cvr_thresh .AND. cvr(k+1) < cvr_thresh) THEN
      nlyr = nlyr + 1
      a(nlyr) = cvr(k)
      ik(nlyr) = k
    ELSE
      IF(nlyr >= 1) THEN
        IF(cvr(k) > a(nlyr)) THEN
          a(nlyr) = cvr(k)         ! Max cldcvr within a layer
          ik(nlyr) = k             ! the lvl index for the max. cldcv
        END IF
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
!  Get temperatures of the maximum cldcvr level in each cld deck.
!
!-----------------------------------------------------------------------
!
  DO n = 1,nlyr
    k = ik(n)
    temp_lyr(n) = t_1d(k)
  END DO ! n

!
!-----------------------------------------------------------------------
!
!  Add a layer for the ground
!
!-----------------------------------------------------------------------
!
  nlyr = nlyr + 1
  a(nlyr) = 1.0
  temp_lyr(nlyr) = t_gnd_k
!
!-----------------------------------------------------------------------
!
!  Convert cloud layer fractions to "cross-section" seen from
!  satellite.  This solves for the "f" array given the "a" array
!
!-----------------------------------------------------------------------
!
  a(1) = MIN(a(1),1.0)
  f(1) = a(1)
  sumf = f(1)
  IF(nlyr >= 2) THEN
    DO n = 2,nlyr
      a(n) = MIN(a(n),1.0)       ! max cldcvr in one cld layer
      f(n) = a(n) * (1.0 - sumf) ! fraction of radiation reaches
                                 ! top of atm from each cld layer
      sumf = sumf + f(n)         ! fraction of radiation blked
                                 ! /attened by all cld lyrs above
    END DO ! n
  END IF ! nlyr
!
!-----------------------------------------------------------------------
!
!  Calculate total radiance from all cloud layers + ground
!
!-----------------------------------------------------------------------
!
  rsum = 0
  DO n = 1,nlyr
    rsum = rsum + temp_to_rad(isat,temp_lyr(n)) * f(n)
  END DO ! n
!
!-----------------------------------------------------------------------
!
!  Convert to effective temperature and compare to
!  the observed brightness temp
!
!-----------------------------------------------------------------------
!
  t_effective = rad_to_temp(isat,rsum)

  RETURN
END SUBROUTINE cvr_to_t11mue


!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE APPLY_CORRECTION               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE apply_correction (cldcv_1dref,cldcv_1d,nz,ilyr,f,delcv)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  To apply the correction to the cloud cover field.
!  The amount for correction is an input parameter.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  05/02/96
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
!
!-----------------------------------------------------------------------
!
!  INPUT:
  INTEGER :: nz        ! number of cloud analysis levels
  REAL :: cldcv_1dref(nz) ! the column cldcvr to be corr.
  INTEGER :: ilyr(nz)  ! cld deck(layer) indices for ea.cldy level
  REAL :: f(nz)           ! fracn of rad. blocked by each cloud deck
  REAL :: delcv            ! the cld cvr increment added to the input
!
!  OUTPUT:
  REAL :: cldcv_1d(nz)     ! corrected cloud cover
!
!  LOCAL:
  INTEGER :: k
!
!-----------------------------------------------------------------------
!
!  Apply correction to 1D cloud cover field
!
!-----------------------------------------------------------------------
!
  DO k = 1,nz
!jz     if( ilyr(k).gt.0 .and. f(ilyr(k)).gt.0.0) then
    IF(ilyr(k) > 0) THEN
      IF(f(ilyr(k)) > 0.0) THEN
        cldcv_1d(k) = MIN(cldcv_1dref(k)+delcv, 1.0)
        cldcv_1d(k) = MAX(cldcv_1d(k), 0.0)
      ELSE
        cldcv_1d(k) = MAX(MIN(cldcv_1dref(k), 1.0), 0.0)
      END IF
    ELSE
      cldcv_1d(k) = MAX(MIN(cldcv_1dref(k), 1.0), 0.0)
    END IF
  END DO

  RETURN
END SUBROUTINE apply_correction
