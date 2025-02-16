!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE get_arps_dims_atts            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_arps_dims_atts(nchr,hinfmt,hisfile,lenfil,               &
                            nx,ny,nz,nzsoil,nstyps,                     &
                            year,month,day,hour,minute,second,time,     &
                            mapproj, sclfct,trulat1,trulat2,trulon,     &
                            ctrlat,ctrlon,dx,dy,ireturn)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Get ARPS dimensions and attributes from ARPS time-independent
!   history file.
!
!-----------------------------------------------------------------------
!
! Author: Yunheng Wang (09/01/2005)
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: nchr
  INTEGER, INTENT(IN)  :: hinfmt
  CHARACTER(*), INTENT(IN)  :: hisfile
  INTEGER, INTENT(IN)  :: lenfil
  INTEGER, INTENT(OUT) :: nx,ny,nz,nzsoil,nstyps
  INTEGER, INTENT(OUT) :: year, month,day,hour,minute,second
  REAL,    INTENT(OUT) :: time
  INTEGER, INTENT(OUT) :: mapproj
  REAL,    INTENT(OUT) :: sclfct
  REAL,    INTENT(OUT) :: trulat1,trulat2,trulon
  REAL,    INTENT(OUT) :: ctrlat,ctrlon
  REAL,    INTENT(OUT) :: dx,dy
  INTEGER, INTENT(OUT) :: ireturn

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: idummy
  REAL    :: rdummy

  CHARACTER(LEN=40) :: fmtverin
  CHARACTER(LEN=10) :: tmunit
  CHARACTER(LEN=80) :: runname
  CHARACTER(LEN=12) :: label
  INTEGER           :: nocmnt, totin
  CHARACTER(LEN=80) :: cmnt(50)
  REAL              :: thisdmp, tstop
  REAL              :: latitud, xgrdorg, ygrdorg, umove, vmove

  REAL, ALLOCATABLE :: x(:)
  REAL, ALLOCATABLE :: y(:)

  INTEGER           :: i

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF( hinfmt == 1 ) THEN

    CALL getunit( nchr )

    OPEN(UNIT=nchr,FILE=hisfile(1:lenfil),STATUS='old',                 &
         FORM='unformatted',IOSTAT=ireturn)

    READ(nchr) fmtverin
    READ(nchr) runname
    READ(nchr) nocmnt
    IF( nocmnt > 0 ) THEN
      DO i=1,nocmnt
        READ(nchr) cmnt(i)
      END DO
    END IF

    READ(nchr) time,tmunit

    READ(nchr) nx, ny, nz,nzsoil

    READ(nchr)      idummy,idummy,idummy,idummy,idummy,                 &
                    idummy,idummy,idummy,idummy,totin,                  &
                    idummy,idummy,idummy,mapproj,month,                 &
                    day, year, hour, minute, second

    READ(nchr)      idummy,idummy,idummy,idummy,idummy,idummy,          &
                    idummy,idummy,idummy,idummy,idummy,idummy,          &
                    idummy,idummy,idummy,idummy,idummy,idummy,          &
                    idummy,idummy,idummy,idummy,idummy,idummy

    READ(nchr)      umove,  vmove, xgrdorg,ygrdorg,trulat1,             &
                    trulat2,trulon,sclfct,rdummy,rdummy,                &
                    rdummy,rdummy,rdummy,rdummy,rdummy,                 &
                    tstop,thisdmp,latitud,ctrlat,ctrlon

    IF ( totin /= 0 ) THEN


      READ(nchr) nstyps,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy

      IF ( nstyps < 1 ) nstyps = 1

      READ(nchr) rdummy,rdummy,rdummy,rdummy,rdummy,                    &
                 rdummy,rdummy,rdummy,rdummy,rdummy,                    &
                 rdummy,rdummy,rdummy,rdummy,rdummy,                    &
                 rdummy,rdummy,rdummy,rdummy,rdummy
    END IF

    ALLOCATE(x(nx), STAT = idummy)
    ALLOCATE(y(ny), STAT = idummy)

    READ(nchr) label
    READ(nchr) x

    READ(nchr) label
    READ(nchr) y

    CLOSE (nchr)
    CALL retunit( nchr )

    dx = x(2) - x(1)
    dy = y(2) - y(1)

    DEALLOCATE(x,y)

  ELSE IF (hinfmt == 3) THEN                  ! HDF 4 format

    CALL hdfopen(hisfile(1:lenfil),1,nchr)

    CALL hdfrdi(nchr,"nx",nx,ireturn)
    CALL hdfrdi(nchr,"ny",ny,ireturn)
    CALL hdfrdi(nchr,"nz",nz,ireturn)
    CALL hdfrdi(nchr,"nzsoil",nzsoil,ireturn)
    CALL hdfrdi(nchr,"nstyp", nstyps,ireturn)

    CALL hdfrdi(nchr,"month", month, ireturn)
    CALL hdfrdi(nchr,"day",   day,   ireturn)
    CALL hdfrdi(nchr,"year",  year,  ireturn)
    CALL hdfrdi(nchr,"hour",  hour,  ireturn)
    CALL hdfrdi(nchr,"minute",minute,ireturn)
    CALL hdfrdi(nchr,"second",second,ireturn)

    CALL hdfrdi(nchr,"mapproj",mapproj,ireturn)
    CALL hdfrdr(nchr,"trulat1",trulat1,ireturn)
    CALL hdfrdr(nchr,"trulat2",trulat2,ireturn)
    CALL hdfrdr(nchr,"trulon", trulon, ireturn)
    CALL hdfrdr(nchr,"sclfct", sclfct, ireturn)
    CALL hdfrdr(nchr,"ctrlat", ctrlat, ireturn)
    CALL hdfrdr(nchr,"ctrlon", ctrlon, ireturn)

    CALL hdfrdr(nchr,"time",time,ireturn)

    CALL hdfrdr(nchr,"dx",dx,ireturn)
    CALL hdfrdr(nchr,"dy",dy,ireturn)

    CALL hdfclose(nchr,ireturn)

  ELSE IF (hinfmt == 7 .OR. hinfmt == 8) THEN ! NetCDF format

    CALL netopen(hisfile(1:lenfil),'R',nchr)

    CALL net_getdims(nchr,nx,ny,nz,nzsoil,nstyps,ireturn)

    CALL netreadatts(nchr,0,runname,nocmnt,cmnt,dx,dy,                  &
                     year,month,day,hour,minute,second,thisdmp,tstop,   &
                     mapproj,sclfct,trulat1,trulat2,trulon,latitud,     &
                     ctrlat,ctrlon,xgrdorg,ygrdorg,umove,vmove,         &
                     rdummy,rdummy,rdummy,rdummy,rdummy,rdummy,         &
                     rdummy,rdummy,rdummy,rdummy,rdummy,                &
                     rdummy,rdummy,rdummy,                              &
                     idummy,idummy,idummy,idummy,idummy,idummy,         &
                     idummy,idummy,idummy,idummy,idummy,                &
                     idummy,idummy,idummy,idummy,                       &
                     idummy,idummy,ireturn)

!    CALL netreadTime(nchr,1,'Time',time)
    time = 0.0

    CALL netclose(nchr)

  ELSE

    WRITE(6,'(a,i3,a)')                                                 &
        ' Data format flag had an invalid value ',                      &
          hinfmt ,' program stopped.'
    CALL arpsstop('arpsstop called from get_arps_dims_atts wrong flag',1)

  END IF

  RETURN
END SUBROUTINE get_arps_dims_atts
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE write_static_attribute        ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_static_attribute(nfid,io_form,start_date,              &
                   grid_id,parent_id,i_parent_start,j_parent_start,     &
                   i_parent_end,j_parent_end,parent_grid_ratio,         &
                   nx,ny,dx,dy,                                         &
                   mapproj,trulat1,trulat2,trulon,ctrlat_moad,          &
                   ctrlat,ctrlon,                                       &
                   lat_ll,lat_ul,lat_ur,lat_lr,                         &
                   lon_ll,lon_ul,lon_ur,lon_lr, istatus)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!   Write attributes of WRF static file.
!
!-----------------------------------------------------------------------
  USE wrf_metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nfid
  INTEGER, INTENT(IN)  :: io_form
  CHARACTER(*), INTENT(IN)  :: start_date
  INTEGER, INTENT(IN)  :: grid_id, parent_id
  INTEGER, INTENT(IN)  :: i_parent_start,j_parent_start
  INTEGER, INTENT(IN)  :: i_parent_end, j_parent_end
  INTEGER, INTENT(IN)  :: parent_grid_ratio
  INTEGER, INTENT(IN)  :: nx,ny
  REAL,    INTENT(IN)  :: dx,dy
  INTEGER, INTENT(IN)  :: mapproj
  REAL,    INTENT(IN)  :: trulat1,trulat2,trulon
  REAL,    INTENT(IN)  :: ctrlat_moad
  REAL,    INTENT(IN)  :: ctrlat, ctrlon
  REAL,    INTENT(IN)  :: lat_ll(4),lat_ul(4),lat_ur(4),lat_lr(4)
  REAL,    INTENT(IN)  :: lon_ll(4),lon_ul(4),lon_ur(4),lon_lr(4)
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER           :: map_proj
  CHARACTER(LEN=20) :: map_string
  REAL              :: latcorns(16), loncorns(16)
  INTEGER           :: i,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  map_string(1:20) = ' '

  IF(ABS(mapproj) == 1) THEN
    map_proj   = 2
    map_string = 'POLAR STEREOGRAPHIC'
  ELSE IF(ABS(mapproj) == 2) THEN
    map_proj   = 1
    map_string = 'LAMBERT CONFORMAL'
  ELSE IF(ABS(mapproj) == 3) THEN
    map_proj = 3
    map_string = 'MERCATOR'
  ELSE
    WRITE(6,*) 'Unknown map projection, ', mapproj
    STOP
  END IF

  DO i = 1,4
    k = (i-1)*4+1
    latcorns(k) = lat_ll(i)
    loncorns(k) = lon_ll(i)
    k = k+1
    latcorns(k) = lat_ul(i)
    loncorns(k) = lon_ul(i)
    k = k+1
    latcorns(k) = lat_ur(i)
    loncorns(k) = lon_ur(i)
    k = k+1
    latcorns(k) = lat_lr(i)
    loncorns(k) = lon_lr(i)
  END DO

  IF (io_form == 7) CALL enter_ncd_define(nfid,istatus)

  CALL put_dom_ti_char(nfid,io_form,'TITLE',                            &
                       'WRF static file from ARPS2WRF',istatus)
  CALL put_dom_ti_char(nfid,io_form,'simulation_name',                  &
                       'WRFSTATIC',istatus)
  CALL put_dom_ti_char(nfid,io_form,'START_DATE',                       &
                       start_date(1:19),istatus)

  CALL put_dom_ti_integer(nfid,io_form,'DYN_OPT',      2,istatus)
  CALL put_dom_ti_integer(nfid,io_form,'WEST-EAST_GRID_DIMENSION',      &
                          nx, istatus)
  CALL put_dom_ti_integer(nfid,io_form,'SOUTH-NORTH_GRID_DIMENSION',    &
                          ny, istatus)

  CALL put_dom_ti_integer(nfid,io_form, 'GRID_ID',       grid_id,   istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'PARENT_ID',     parent_id, istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'I_PARENT_START',i_parent_start, istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'J_PARENT_START',j_parent_start, istatus)

  CALL put_dom_ti_integer(nfid,io_form, 'I_PARENT_END',  i_parent_end,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'J_PARENT_END',  i_parent_end,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'PARENT_GRID_RATIO',parent_grid_ratio,istatus)

  CALL put_dom_ti_real   (nfid,io_form, 'DX',dx,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'DY',dy,istatus)

  CALL put_dom_ti_varreal(nfid,io_form, 'corner_lats',latcorns,16,istatus)
  CALL put_dom_ti_varreal(nfid,io_form, 'corner_lons',loncorns,16,istatus)

  CALL put_dom_ti_real   (nfid,io_form, 'CEN_LAT',ctrlat,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'CEN_LON',ctrlon,istatus)

  CALL put_dom_ti_integer(nfid,io_form, 'FLAG_STATIC',1,istatus)

  CALL put_dom_ti_char(nfid,io_form,'map_projection',TRIM(map_string),  &
                       istatus)

  CALL put_dom_ti_integer(nfid,io_form, 'MAP_PROJ',     map_proj,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'MOAD_CEN_LAT', ctrlat_moad,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'STAND_LON',    trulon,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'TRUELAT1',     trulat1,istatus)
  CALL put_dom_ti_real   (nfid,io_form, 'TRUELAT2',     trulat2,istatus)

  CALL put_dom_ti_integer(nfid,io_form, 'ISWATER',ISWATER,istatus)
  CALL put_dom_ti_integer(nfid,io_form, 'ISICE',  ISICE,  istatus)
  CALL put_dom_ti_char   (nfid,io_form, 'MMINLU', 'USGS', istatus)

  IF (io_form == 7) CALL exit_ncd_define(nfid,istatus)

  RETURN
END SUBROUTINE write_static_attribute
!
! ================================================================
!
SUBROUTINE RDGEODAT(n2,n3,xt,yt,deltax,deltay,                          &
                  ofn,wvln,silwt,maxdatacat,                            &
                  datr,dats,datln,datlt,which_data, istat)

  IMPLICIT NONE

  INTEGER, INTENT(IN) ::  n2,n3
  REAL,    INTENT(IN) ::  deltax,deltay,wvln,silwt
              ! wvln = TOPTWVL_PARM_WRF from wrfsi.nl section &hgridspec
              ! siwt = SILAVWT_PARM_WRF
  REAL,    INTENT(IN) ::  xt(n2),yt(n3)
  INTEGER, INTENT(IN) ::  maxdatacat
  CHARACTER(LEN=*),INTENT(IN) ::  OFN

  REAL,    INTENT(OUT):: datr(n2,n3)
  REAL,    INTENT(OUT):: dats(n2,n3,maxdatacat)
  REAL,    INTENT(OUT):: datln(n2,n3)
  REAL,    INTENT(OUT):: datlt(n2,n3)

  LOGICAL, INTENT(OUT)::  which_data

  INTEGER, INTENT(OUT):: istat

!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: IODIM = 5800000
    ! this parameter can be increased if we ever read data with
    ! resolution finer than 30 sec, or if the tilesize for 30s
    ! data becomes greater than 10x10 deg.
    !
  INTEGER :: lb,lbf
  INTEGER :: mof,np,niq,njq,nx,ny

  REAL    :: deltallo,deltaxq,deltayq, deltaxp,deltayp
  INTEGER :: no, iblksizo, isbego, iwbego
  REAL    :: rsoff,rwoff

  INTEGER :: lcat

  CHARACTER(LEN=180) ::TITLE

  REAL            :: std_lon = -100.0    ! not used anywhere
  REAL, PARAMETER :: erad    = 6371200.0 ! not used anywhere

  CHARACTER(LEN=180) :: tmpstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! *********************
  nx = n2-1
  ny = n3-1
! **********************

  lcat  = 1
  lbf   = LEN_TRIM(ofn)

  tmpstr(1:180) = ' '
  write(tmpstr,'(a)') ofn(1:lbf)

  title = ofn(1:lbf)//'HEADER'
  lb = LEN_TRIM(title)

  OPEN(29,STATUS='OLD',FILE=title(1:lb),FORM='FORMATTED')

  READ(29,*) iblksizo,no,isbego,iwbego,rsoff,rwoff
  print *,'title         = ',title
  print *,'rsoff,rwoff   = ',rsoff,rwoff
  print *,'isbego,iwbego = ',isbego,iwbego
  print *,'iblksizo,no   = ',iblksizo,no
  CLOSE(29)

  IF (no > 1) THEN
     deltallo = FLOAT(iblksizo)/FLOAT(no-1)
  ELSE IF(no == 1) THEN
     deltallo = FLOAT(iblksizo)/FLOAT(no)
  ELSE
     print*,'HEADER value NO = 0'
     istat = -1
     RETURN
  ENDIF

  mof = iodim/(no*no)

  IF (ofn(lbf:lbf) == 'G' .OR. ofn(lbf:lbf) == 'A') THEN  ! Green_frac or Albedo
    lcat=12
    IF (no == 1250) mof = 1
  END IF

! MOF determines the number of files held in buffer while reading
! DEM data; it saves some time when buffer data can be used instead
! of reading DEM file again. Originally MOF was 4.

  IF (mof > 10) mof = 5

  deltaxq = 0.5*wvln*deltax
  deltayq = 0.5*wvln*deltay
  print *,'deltaxq,deltayq=',deltaxq,deltayq

  niq = INT(FLOAT(nx)*deltax/deltaxq) + 4
  njq = INT(FLOAT(ny)*deltay/deltayq) + 4

  np  = MIN(10,MAX(1,INT(deltaxq/(deltallo*111000.))))
  print *,' np=',np

  deltaxp = deltaxq/FLOAT(np)
  deltayp = deltayq/FLOAT(np)

  CALL SFCOPQR(NO,MOF,NP,NIQ,NJQ,N2,N3,lcat,                            &
         XT,YT,90.,std_lon,ERAD,RWOFF,RSOFF,                            &
         DELTALLO,DELTAXP,DELTAYP,DELTAXQ,DELTAYQ,                      &
         IBLKSIZO,ISBEGO,IWBEGO,DATR,DATS,DATLN,DATLT,                  &
         tmpstr,WVLN,SILWT,which_data,maxdatacat,istat)

  RETURN
END SUBROUTINE rdgeodat

SUBROUTINE get_path_to_soiltemp_1deg (path_soiltemp_1deg,istatus)

  IMPLICIT NONE
  CHARACTER*256 :: path_soiltemp_1deg
  INTEGER       :: istatus

  CHARACTER(LEN=256) :: path_to_soiltemp_1deg
  COMMON /lapscmn/ path_to_soiltemp_1deg

  path_soiltemp_1deg = path_to_soiltemp_1deg

  RETURN
END SUBROUTINE

SUBROUTINE interpolate_masked_val(nx_src, ny_src, lmask_src, data_src, &
                                 nx_out, ny_out, lmask_out, data_out, &
                                 isrcr, jsrcr, make_srcmask, &
                                 min_val, max_val, def_val, val_mask, &
                                 method)

    IMPLICIT none

    INTEGER, INTENT(IN)           :: nx_src
    INTEGER, INTENT(IN)           :: ny_src
    REAL, INTENT(INOUT)           :: lmask_src(nx_src,ny_src)
    REAL, INTENT(IN)              :: data_src(nx_src,ny_src)
    INTEGER, INTENT(IN)           :: nx_out, ny_out
    REAL, INTENT(IN)              :: lmask_out(nx_out,ny_out)
    REAL, INTENT(INOUT)             :: data_out(nx_out,ny_out)
    REAL, INTENT(IN)              :: isrcr(nx_out,ny_out)
    REAL, INTENT(IN)              :: jsrcr(nx_out,ny_out)
    LOGICAL, INTENT(IN)           :: make_srcmask
    REAL, INTENT(IN)              :: min_val, max_val, def_val, val_mask
    INTEGER, INTENT(IN)           :: method

    INTEGER                       :: i,j,k,l,ilo,jlo, ii, jj, kk, ll
    INTEGER                       :: search_rad
    INTEGER                       :: search_rad_max
    REAL, PARAMETER               :: bad_point = -99999.
    LOGICAL                       :: bad_points_flag
    LOGICAL                       :: fixed_bad
    INTEGER                       :: total_points
    INTEGER                       :: points_val_mask
    INTEGER                       :: points_skipped
    INTEGER                       :: num_to_fix
    INTEGER                       :: fixed_with_out
    INTEGER                       :: fixed_with_src
    INTEGER                       :: fixed_with_def
    INTEGER                       :: close_pts
    INTEGER                       :: ic,jc,iic,jjc
    REAL                          :: deltagx, deltagy
    REAL                          :: dix, djy
    REAL                          :: data_close(4)
    REAL                          :: distance(4)
    REAL                          :: distx,disty
    REAL                          :: sum_dist
    REAL                          :: a,b,c,d,e,f,g,h
    REAL                          :: stl(4,4)
    REAL                          :: valb
    REAL                          :: srcmask(nx_src,ny_src)

    LOGICAL                       :: wrapped_source

  INTEGER, PARAMETER          :: METHOD_NEAREST = 0
  INTEGER, PARAMETER          :: METHOD_LINEAR  = 1
  INTEGER, PARAMETER          :: METHOD_HIGHER = 2

    REAL                      :: oned

    IF ((MINVAL(isrcr).LT.1.0).OR. (MAXVAL(isrcr).GT.nx_src))THEN
      wrapped_source = .TRUE.
    ELSE
      wrapped_source = .FALSE.
    ENDIF
    ! Can we use the source data land mask that was input, or do we need to make it
    ! from the min_val/max_val arguments?

    IF (make_srcmask) THEN
      PRINT '(A)', 'INTERPOLATE_MASKED_VAL: Making source data land mask...'

      ! This code has the flexibility to handled either water data points (lmask =0)
      ! or land points (lmask = 1). So if we are making the landmask using valid
      ! range, but we do not know anything about the variable other than the
      ! valid mask value, we need to figure out which points are 0 and which should be
      ! 1.  To do this, initialize the array to the invalid value, which we determine
      ! from the valid value.
      IF (val_mask .EQ. 0) THEN
        lmask_src(:,:) = 1
      ELSEIF (val_mask .EQ. 1) THEN
        lmask_src(:,:) = 0
      ELSE
        lmask_src(:,:) = -1
      ENDIF

      ! Now figure out where to set the mask using a WHERE statement
      ! NOTE:  In this next line, if val_mask is 2, then the lmask_src
      ! is going to be set to 2, so we need to be careful in some of the
      ! subsequent IF statements where lmask_src is compared to lmask_out
      WHERE((data_src .GE. min_val).AND.(data_src .LE. max_val)) lmask_src = val_mask
    ELSE
      PRINT '(A)', 'INTERPOLATE_MASKED: Using source landmask field.'
    ENDIF


    bad_points_flag = .false.

    ! Initialize counters
    total_points = 0
    points_val_mask = 0
    points_skipped = 0
    num_to_fix = 0
    fixed_with_out = 0
    fixed_with_src = 0
    fixed_with_def = 0

    ! Select interpolation method.  Putting the case statement out here
    ! increases the amount of replicated code but should be more efficient
    ! than checking this condition at every grid point.

    SELECT CASE(method)

      CASE(METHOD_NEAREST)
        ! Use nearest neigbor
        PRINT '(A)', 'Masked interpolation using nearest neighbor value...'
        out_j_loop_1: DO j = 1, ny_out
          out_i_loop_1: DO i = 1, nx_out

            total_points = total_points + 1
            ! We only need to process this point if the lmask_out is equal
            ! to val_mask.  For example, if one is processing soil parameters,
            ! which are only valid at land points (lmask = 1), then the user
            ! passes in 1. for the val_mask.  Any output point that is water
            ! (lmask = 0) will then be skipped.  During this loop, if we have
            ! points that cannot be properly assigned, we will mark them as bad
            ! and set the bad_points_lag.  Exception is if val_mask = 2, which
            ! implies that we are doing masked interpolation for both land
            ! and water, but allowing only water points to influence water
            ! points and land points to influence landpoints.

            IF ((lmask_out(i,j) .EQ. val_mask).OR.(val_mask .EQ. 2)) THEN
              ! Process this point
              points_val_mask = points_val_mask + 1
              ilo = NINT(isrcr(i,j))
              ! Account for horizontal wrapping as in the case
              ! of global data sets.  This is the only time
              ! it is possible for ilo < 1 or ilo > nx_src
              IF (ilo .LT.1) ilo = nx_src
              IF (ilo .GT. nx_src) ilo = 1

              jlo = NINT(jsrcr(i,j))

              ! See if this point can be used
              IF ((lmask_src(ilo,jlo).EQ. lmask_out(i,j)).OR. &
                  (lmask_src(ilo,jlo).EQ.2)) THEN
                data_out(i,j) = data_src(ilo,jlo)
              ELSE
                data_out(i,j) = bad_point
                num_to_fix = num_to_fix + 1
                bad_points_flag = .true.
              ENDIF
            ELSE
              ! The output grid does not require a value for this point
              ! But do not zero out in case this is a field begin
              ! done twice (once for water and once for land, e.g.
              ! SKINTEMP
              !data_out(i,j) = 0.
              points_skipped = points_skipped + 1
            ENDIF
          ENDDO out_i_loop_1
        ENDDO out_j_loop_1

      CASE (METHOD_LINEAR)
        ! Use a 4-point interpolation
        PRINT '(A)', 'Masked interpolation using 4-pt linear interpolation...'
        deltagx = 1.
        deltagy = 1.
        out_j_loop_2: DO j = 1, ny_out
          out_i_loop_2: DO i = 1, nx_out

            total_points = total_points + 1
            ! We only need to process this point if the lmask_out is equal
            ! to val_mask.  For example, if one is processing soil parameters,
            ! which are only valid at land points (lmask = 1), then the user
            ! passes in 1. for the val_mask.  Any output point that is water
            ! (lmask = 0) will then be skipped.  During this loop, if we have
            ! points that cannot be properly assigned, we will mark them as bad
            ! and set the bad_points_lag.

            IF ((lmask_out(i,j) .EQ. val_mask).OR.(val_mask .EQ. 2)) THEN
              ! Process this point
              points_val_mask = points_val_mask + 1

              ! If ilo < 1 or > nx_src, this is a wrapped source dataset

              ilo = FLOOR(isrcr(i,j))
              IF (ilo .EQ. 0) ilo = nx_src
              jlo = MIN(FLOOR(jsrcr(i,j)),ny_src-1)
              dix = isrcr(i,j) - FLOAT(ilo)
              IF (dix .LT.0) dix = dix + FLOAT(nx_src)
              djy = jsrcr(i,j) - FLOAT(jlo)

              ! Loop around the four surrounding points
              ! and count up the number of points we can use based
              ! on common mask value
              close_pts = 0
              sum_dist = 0.
              outer_four_j: DO jc = 0,1
                outer_four_i: DO ic = 0,1
                  iic = ilo + ic
                  IF(iic .GT. nx_src) iic = iic - nx_src
                  jjc = jlo + jc
                  IF ((lmask_src(iic,jjc).EQ. lmask_out(i,j)).OR. &
                      (lmask_src(iic,jjc).EQ.2)) THEN
                    close_pts = close_pts + 1
                    data_close(close_pts) = data_src(iic,jjc)

                    ! Compute distance to this valid point
                    ! in grid units and add to sum of distances (actually,
                    ! we are doing a weight, which is inversely proportional
                    ! to distance)
                    IF (ic .EQ. 0) THEN
                      distx = deltagx - dix
                    ELSE
                      distx =  dix
                    ENDIF
                    IF (jc .EQ. 0) THEN
                      disty = deltagy - djy
                    ELSE
                      disty = djy
                    ENDIF
                    distance(close_pts) = SQRT(distx**2+disty**2)
                    sum_dist = sum_dist + distance(close_pts)
                  ENDIF
                ENDDO outer_four_i
              ENDDO outer_four_j

              ! Did we find at least one point in the surrounding four
              ! that was usable?

              IF (close_pts .GT. 0) THEN

                ! If we have all four points, then do bilinear interpolation
                IF (close_pts .EQ. 4) THEN
                   data_out(i,j) = ((deltagx - dix)*((deltagy-djy)*data_close(1) &
                                 + djy*data_close(3)) &
                                 + dix*((deltagy-djy)*data_close(2) &
                                 + djy*data_close(4))) &
                                 / (deltagx * deltagy)
                ELSE IF ((close_pts .GT. 1).AND.(close_pts .LT. 4)) THEN

                  ! Simple distance-weighted average by computing
                  ! the sum of all distances to each point and using
                  ! each individual distance divided by the total
                  ! distance as the weighting

                  data_out(i,j) = 0.
                  DO k = 1, close_pts
                    data_out(i,j) = data_out(i,j) + &
                                    (distance(k)/sum_dist) * data_close(k)
                  ENDDO
                ELSE
                  ! Set output value = to one point we found
                  data_out(i,j) = data_close(1)
                ENDIF
              ELSE
                bad_points_flag = .true.
                data_out(i,j) = bad_point
                num_to_fix = num_to_fix + 1
              ENDIF
            ELSE
              ! The output grid does not require a value for this point
              ! But do not zero out in case this is a field begin
              ! done twice (once for water and once for land, e.g.
              ! SKINTEMP
              points_skipped = points_skipped + 1
            ENDIF
          ENDDO out_i_loop_2
        ENDDO out_j_loop_2

      CASE (METHOD_HIGHER)
        ! 16-point interpolation
        out_j_loop_3: DO j = 1, ny_out
          out_i_loop_3: DO i = 1, nx_out

            total_points = total_points + 1
            ! We only need to process this point if the lmask_out is equal
            ! to val_mask.  For example, if one is processing soil parameters,
            ! which are only valid at land points (lmask = 1), then the user
            ! passes in 1. for the val_mask.  Any output point that is water
            ! (lmask = 0) will then be skipped.  During this loop, if we have
            ! points that cannot be properly assigned, we will mark them as bad
            ! and set the bad_points_lag.

            IF ((lmask_out(i,j) .EQ. val_mask).OR.(val_mask .EQ. 2)) THEN
              ! Process this point
              points_val_mask = points_val_mask + 1

              ! Do a 4x4 loop around in the input data around the output
              ! point to get our 16 points of influence.  Only use those
              ! that are appropriately masked
              valb = 0.
              close_pts = 0
              ilo = INT(isrcr(i,j)+0.00001)
              jlo = INT(jsrcr(i,j)+0.00001)
              dix = isrcr(i,j) - ilo
              djy = jsrcr(i,j) - jlo
              IF ( (ABS(dix).GT.0.0001).OR.(ABS(djy).GT.0.0001) ) THEN
                ! Do the interpolation loop
                stl(:,:) = 0.
                loop_16_1: DO k = 1,4
                  kk = ilo + k - 2
                  IF ((kk .LT. 1).OR.(kk .GT. nx_src)) THEN
                    IF (.NOT. wrapped_source) THEN
                      CYCLE loop_16_1
                    ELSE
                      IF (kk .LT. 1) THEN
                        kk = kk + nx_src
                      ELSE
                        kk = kk - nx_src
                      ENDIF
                    ENDIF
                  ENDIF
                  loop_16_2: DO l = 1, 4
                    ll = jlo + l - 2
                    IF ((ll .LT. 1).OR.(ll .GT. ny_src)) CYCLE loop_16_2

                    ! Check land mask at this source point
                    IF ((lmask_src(kk,ll).NE. lmask_out(i,j)).AND. &
                        (lmask_src(kk,ll).NE.2))  CYCLE loop_16_2
                    ! If we are here, then mask tests passed
                    stl(k,l) = data_src(kk,ll)
                    IF ( (stl(k,l) .EQ. 0.).AND.(min_val.LE.0.).AND. &
                                                (max_val.GE.0.) ) THEN
                      stl = 1.E-5
                    ENDIF
                    close_pts = close_pts + 1
                  ENDDO loop_16_2
                ENDDO loop_16_1

                ! Did we find any valid points?

                IF ( (close_pts .GT. 0).AND. ( &
                  (stl(2,2).GT.0.).AND.(stl(2,3).GT.0.).AND. &
                  (stl(3,2).GT.0.).AND.(stl(3,3).GT.0.)  ) ) THEN
                  a = oned(dix,stl(1,1),stl(2,1),stl(3,1),stl(4,1))
                  b = oned(dix,stl(1,2),stl(2,2),stl(3,2),stl(4,2))
                  c = oned(dix,stl(1,3),stl(2,3),stl(3,3),stl(4,3))
                  d = oned(dix,stl(1,4),stl(2,4),stl(3,4),stl(4,4))
                  valb = oned(djy,a,b,c,d)
                  IF (close_pts .NE. 16) THEN
                    e = oned(djy,stl(1,1),stl(1,2),stl(1,3),stl(1,4))
                    f = oned(djy,stl(2,1),stl(2,2),stl(2,3),stl(2,4))
                    g = oned(djy,stl(3,1),stl(3,2),stl(3,3),stl(3,4))
                    h = oned(djy,stl(4,1),stl(4,2),stl(4,3),stl(4,4))
                    valb = (valb+oned(dix,e,f,g,h)) * 0.5
                  ENDIF
                  data_out(i,j) = valb

                ELSE
                  bad_points_flag = .true.
                  data_out(i,j) = bad_point
                  num_to_fix = num_to_fix + 1
                ENDIF
              ELSE
                ! We are right on a source point, so try to use it
                IF ((lmask_src(ilo,jlo).EQ.lmask_out(i,j)).OR.&
                    (lmask_src(ilo,jlo).EQ.2)) THEN
                  data_out(i,j) = data_src(ilo,jlo)
                ELSE
                  bad_points_flag = .true.
                  data_out(i,j) = bad_point
                  num_to_fix = num_to_fix + 1
                ENDIF
              ENDIF
            ELSE
              ! The output grid does not require a value for this point
              !data_out(i,j) = 0.
              points_skipped = points_skipped + 1
            ENDIF
          ENDDO out_i_loop_3
        ENDDO out_j_loop_3

    END SELECT

     ! Do we need to correct bad points?

    IF (bad_points_flag) THEN

      search_rad_max = 10
      fix_bad_j: DO j = 1, ny_out
        fix_bad_i: DO i = 1, nx_out

          IF (data_out(i,j).NE. bad_point) CYCLE fix_bad_i

          ! First, search for nearest non-bad point in the output domain
          ! which is usually higher resolution.
          fixed_bad = .false.
          search_out_loop: DO search_rad = 1, search_rad_max
            search_out_j: DO ll = -(search_rad-1), (search_rad-1),1
              jj = j + ll
              IF ((jj .LT. 1).OR.(jj .GT. ny_out)) CYCLE search_out_j
              search_out_i: DO kk = -(search_rad), search_rad, 1
                 ii = i + kk
                 IF ((ii .LT. 1).OR.(ii .GT. nx_out)) CYCLE search_out_i
                 IF ((data_out(ii,jj).NE.bad_point).AND. &
                    (lmask_out(ii,jj) .EQ. lmask_out(i,j)) ) THEN
                  data_out(i,j) = data_out(ii,jj)
                  fixed_bad = .true.
                  fixed_with_out = fixed_with_out + 1
                  EXIT search_out_loop
                ENDIF
              ENDDO search_out_i
            ENDDO search_out_j
          ENDDO search_out_loop

          ! Did we fix the point?  If not, then do same search on src data.
          IF (.NOT. fixed_bad) THEN
            search_rad_max = 10
            search_src_loop: DO search_rad = 1, search_rad_max
              search_src_j: DO ll = -(search_rad-1), (search_rad-1),1
                jj = NINT(jsrcr(i,j)) + ll
                IF ((jj .LT. 1).OR.(jj .GT. ny_src)) CYCLE search_src_j
                search_src_i: DO kk = -(search_rad), search_rad, 1
                   ii = NINT(isrcr(i,j)) + kk
                   IF ((ii .LT. 1).OR.(ii .GT. nx_src)) THEN
                     IF (.NOT. wrapped_source) THEN
                       CYCLE search_src_i
                     ELSE
                       IF (ii.LT.1) THEN
                         ii = nx_src + ii
                       ELSE
                         ii = ii - nx_src
                       ENDIF
                     ENDIF
                   ENDIF
                   IF ((lmask_src(ii,jj).EQ.lmask_out(i,j)).OR. &
                       (lmask_src(ii,jj).EQ.2)) THEN
                     data_out(i,j) = data_src(ii,jj)
                     fixed_bad = .true.
                     fixed_with_src = fixed_with_src + 1
                     EXIT search_src_loop
                   ENDIF
                ENDDO search_src_i
              ENDDO search_src_j
            ENDDO search_src_loop
          ENDIF
          ! Now is the point fixed?  If not, we have to use a default value.
          IF (.NOT.fixed_bad) THEN
            fixed_with_def = fixed_with_def + 1
            data_out(i,j) = def_val
            PRINT '(A,F10.3,A,2I5)', 'INTERPOLATE_MASKED: Bogus value of ', def_val, &
                ' used at point ', i, j
          ENDIF

        ENDDO fix_bad_i
      ENDDO fix_bad_j
    ENDIF
    PRINT '(A)',     '----------------------------------------'
    PRINT '(A)',     'MASKED INTERPOLATION SUMMARY: '
    PRINT '(A,I10)', '  TOTAL POINTS IN GRID:       ', total_points
    PRINT '(A,I10)', '  POINTS NEEDING VALUES:      ', points_val_mask
    PRINT '(A,I10)', '  POINTS NOT REQUIRED:        ', points_skipped
    PRINT '(A,I10)', '  POINTS NEEDING FIX:         ', num_to_fix
    PRINT '(A,I10)', '  POINTS FIXED WITH OUT GRID: ', fixed_with_out
    PRINT '(A,I10)', '  POINTS FIXED WITH SRC GRID: ', fixed_with_src
    PRINT '(A,I10)', '  POINTS FIXED WITH DEF VAL:  ', fixed_with_def
    PRINT '(A)',     '----------------------------------------'
    RETURN
END SUBROUTINE interpolate_masked_val

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION oned(x,a,b,c,d)

      IMPLICIT NONE

      REAL :: x,a,b,c,d,oned

      oned = 0.

      IF      ( x .EQ. 0. ) THEN
         oned = b
      ELSE IF ( x .EQ. 1. ) THEN
         oned = c
      END IF

      IF(b*c.NE.0.) THEN
         IF ( a*d .EQ. 0. ) THEN
            IF      ( ( a .EQ. 0 ) .AND. ( d .EQ. 0 ) ) THEN
               oned = b*(1.0-x)+c*x
            ELSE IF ( a .NE. 0. ) THEN
               oned = b+x*(0.5*(c-a)+x*(0.5*(c+a)-b))
            ELSE IF ( d .NE. 0. ) THEN
               oned = c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c))
            END IF
         ELSE
            oned = (1.0-x)*(b+x*(0.5*(c-a)+x*(0.5*(c+a)-b)))+ &
                   x*(c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c)))
         END IF
      END IF

END FUNCTION oned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE latlon_2_llij(np,glat,glon,lli,llj,lat0,lon0,dlat,dlon,cgrddef)

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: np
  REAL,    INTENT(IN)  :: glat(np), glon(np)
  REAL,    INTENT(OUT) :: lli(np),  llj(np)
  REAL,    INTENT(IN)  :: lat0, lon0
  REAL,    INTENT(IN)  :: dlat, dlon
  CHARACTER(LEN=1), INTENT(IN) :: cgrddef

  INTEGER :: n
  REAL    :: diff
  REAL    :: dlond, dlatd

  dlond=dlon
  dlatd=dlat
  IF (cgrddef .EQ. 'S') THEN
    DO n=1,np
      diff=glon(n)-lon0
      IF (diff <    0.) diff=diff+360.
      IF (diff >= 360.) diff=diff-360.
      lli(n) = diff/dlond+1.
      llj(n) = (glat(n)-lat0)/dlatd+1.
    END DO

  ELSE IF (cgrddef .EQ. 'N') THEN
    DO n=1,np
      diff=glon(n)-lon0
      IF (diff <    0.) diff=diff+360.
      IF (diff >= 360.) diff=diff-360.
      lli(n) = diff/dlon+1.
      llj(n) = (lat0-glat(n))/dlat+1.
    END DO

  ELSE
    print*, 'you must specify whether the standard  &
           & lat is Southern or Northern boundary'
  END IF

  RETURN
END SUBROUTINE latlon_2_llij

FUNCTION bilinear_interp(i,j,imax,jmax,array_2d) RESULT (zo)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j,imax,jmax
  REAL,    INTENT(IN) :: array_2d(imax,jmax)
  REAL                :: zo

  REAL ::  z1,z2,z3,z4
  REAL ::  fraci,fracj

  fraci = 0.5
  fracj = 0.5

  z1 = array_2d(i  , j  )
  z2 = array_2d(i-1, j  )
  z3 = array_2d(i-1, j-1)
  z4 = array_2d(i  , j-1)

  zo = Z1 + (Z2-Z1)*fraci + (Z4-Z1)*fracj - (Z2+Z4-Z3-Z1)*fraci*fracj

END FUNCTION bilinear_interp
