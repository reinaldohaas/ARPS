!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTSOIL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtsoil(nx,ny,nzsoil,nstyps,soiloutfl,dx,dy,zpsoil,          &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   zpsoilout,tsoilout,qsoilout,wcanpout,snowdout,       &
                   tsoil,qsoil,wetcanp,snowdpth,soiltyp )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write the soil model variables into file soiloutfl.
!
!  NOTE:
!
!    Any changes made in this subroutine should also be made in
!    wrtjoinsoil below.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  04/20/1995
!
!  MODIFICATION HISTORY:
!
!  08/07/1995 (M. Xue)
!  Added file name to the argument list.
!
!  2000/03/29 (Gene Bassett)
!  Removed the globcst.inc include.
!
!  2000/03/29 (Gene Bassett)
!  Added HDF4 format.
!
!  2000/10/27 (Gene Bassett)
!  Added soiltyp to soil file to track consistency of soil data.
!
!  05/13/2002 (J. Brotzge)
!  Modified arrays, call statements to allow for multiple soil schemes.
!
!  08/16/2004 (Yunheng Wang)
!  Added NetCDF format.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx        Number of model grid points in the x-dir. (east/west)
!  ny        Number of model grid points in the y-dir. (north/south)
!  nzsoil    Number of model grid points in the soil.
!  nstyps    Number of soil types.
!
!  soiloutfl Name of the soil model data output file
!
!  zpsoilout Flag for output of soil level depth
!  tsoilout  Flag for output of soil temperature
!  qsoilout  Flag for output of soil moisture
!  wcanpout  Flag for output of Canopy water amount
!
!  tsoil     Soil temperature (K)
!  qsoil     Soil moisture (m3/m3)
!  wetcanp   Canopy water amount
!  snowdpth  Snow depth (m)
!  soiltyp   Soil type in model domain
!
!  OUTPUT:
!
!  Output written to the surface data file, soilinfl:
!
!  nx       Number of model grid points in the x-direction
!  ny       Number of model grid points in the y-direction
!  nzsoil   Number of model grid points in the soil
!
!  mapproj  Type of map projection used to setup the analysis grid.
!  trulat1  1st real true latitude of map projection.
!  trulat2  2nd real true latitude of map projection.
!  trulon   Real true longitude of map projection.
!  sclfct   Map scale factor. At latitude = trulat1 and trulat2
!
!  dx       Model grid spacing in the x-direction east-west (meters)
!  dy       Model grid spacing in the y-direction east-west (meters)
!  ctrlat   Lat. at the origin of the model grid (deg. N)
!  ctrlon   Lon. at the origin of the model grid (deg. E)
!
!  tsoilout Flag for output of soil temperature
!  qsoilout Flag for output of soil moisture
!  wcanpout Flag for output of Canopy water amount
!
!  zpsoil   Soil level depth (m)
!  tsoil    Soil temperature (K)
!  qsoil    Soil moisture    (m3/m3)
!  wetcanp  Canopy water amount
!  snowdpth Snow depth       (m)
!  soiltyp  Soil type in model domain
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx            ! Number of grid points in the x-direction
  INTEGER :: ny            ! Number of grid points in the y-direction
  INTEGER :: nzsoil        ! Number of grid points in the soil

  INTEGER :: nstyps        ! Number of soil types at each grid point
  CHARACTER(LEN=*) :: soiloutfl ! Name of the output file

  REAL :: dx
  REAL :: dy
  INTEGER :: mapproj       ! Map projection
  REAL :: trulat1          ! 1st real true latitude of map projection
  REAL :: trulat2          ! 2nd real true latitude of map projection
  REAL :: trulon           ! Real true longitude of map projection
  REAL :: sclfct           ! Map scale factor
  REAL :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon           ! Center longitude of the model domain (deg. E)

  INTEGER :: zpsoilout     ! Flag for output of zpsoil
  INTEGER :: tsoilout      ! Flag for output of tsoil
  INTEGER :: qsoilout      ! Flag for output of qsoil
  INTEGER :: wcanpout      ! Flag for output of wetcanp
  INTEGER :: stypout       ! Flag for output of soiltyp (set to 1 if any of
                           ! the above flags are set)
  INTEGER :: snowdout      ! Flag for output of snowdpth

  REAL :: zpsoil (nx,ny,nzsoil)          ! Soil depth (m)
  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m3/m3)
  REAL :: wetcanp(nx,ny,0:nstyps)        ! Canopy water amount

  REAL :: snowdpth(nx,ny)            ! Snow depth (m)

  INTEGER :: soiltyp (nx,ny,nstyps)  ! Soil type in model domain

  REAL,    ALLOCATABLE :: var3d (:,:,:)
  INTEGER, ALLOCATABLE :: var3di(:,:,:)
  REAL,    ALLOCATABLE :: var4d (:,:,:,:)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,is
  INTEGER :: flunit
  INTEGER :: idummy
  REAL    :: rdummy
  INTEGER :: ierr

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
  REAL :: atmp1(1),atmp2(1)
         ! unused arrays in hdf routines since NO COMPRESSION

  INTEGER :: istat, sd_id

!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver500
  INTEGER  :: intver,intver500

  PARAMETER (fmtver500='* 005.00 GrADS Soilvar Data',intver500=500)

  INTEGER :: npxout, npyout
  INTEGER :: nxout,  nyout
  INTEGER :: ipx,    jpy
  INTEGER :: ia,     ja

  INTEGER :: lenbase
  CHARACTER(LEN=256) :: outfilename

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  idummy = 0
  rdummy = 0.0

  IF (soildmp == 0) RETURN

  IF ( tsoilout+qsoilout+wcanpout /= 0 ) THEN
    stypout = 1
  ELSE
    stypout = 0
  END IF

  IF ( (nproc_x_out > nproc_x .OR. nproc_y_out > nproc_y) .AND.         &
       splitsoil > 0  ) THEN

    npxout = nproc_x_out/nproc_x
    npyout = nproc_y_out/nproc_y

    nxout = (nx-3)/npxout + 3
    nyout = (ny-3)/npyout + 3

  ELSE
    npxout = 1
    npyout = 1
    nxout  = nx
    nyout  = ny
  END IF

  !lenbase = INDEX(soiloutfl,'_',.TRUE.)
  !idummy  = INDEX(soiloutfl,'.',.TRUE.)
  !IF (lenbase < idummy)
  lenbase = LEN_TRIM(soiloutfl)

!
! To avoid nstyp > nstyps, in such a case tsoil, qsoil, wetcanp will be
! undefined for soil types diemension great than nstyps. -- WYH.
!
  nstyp = MIN(nstyp, nstyps)

  IF (myproc == 0)  &
    WRITE (6,'(1x,/,2a)') 'Write soil initial data to ',TRIM(soiloutfl)

!-----------------------------------------------------------------------
!
!  Write out in Fortran unformatted.
!
!-----------------------------------------------------------------------

  IF (soildmp == 1) THEN

    CALL getunit( flunit )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(soiloutfl, '-F f77 -N ieee', ierr)

    intver = intver500  !  for the time being, in the future, we will
                        !  allow to dump data in the different version
                        !  intver will be assigned from input file
    IF (intver == intver500) THEN
      fmtver=fmtver500
    ELSE
      WRITE(6,'(/1x,a,i10,a/)')                        &
        'Data format, intver=',intver,', not found. The Job stopped.'
      CALL arpsstop('arpstop called from wrtsoil.',1)
    END IF

    CALL fnversn(soiloutfl, lenbase)

    OPEN (UNIT=flunit, FILE=trim(soiloutfl),STATUS='unknown',        &
          FORM='unformatted', ACCESS='sequential')

    WRITE (flunit) fmtver

    WRITE (flunit) nx, ny, nzsoil

    WRITE (flunit) mapproj,tsoilout, qsoilout,                       &
                   wcanpout,      0,snowdout,  stypout, zpsoilout,   &
                   idummy,   idummy,  idummy,   idummy,  idummy,     &
                   idummy,   idummy,  idummy,   idummy,  nstyp

    WRITE (flunit)    dx,     dy, ctrlat, ctrlon,trulat1,            &
                 trulat2, trulon, sclfct, rdummy, rdummy,            &
                 rdummy,  rdummy, rdummy, rdummy, rdummy,            &
                 rdummy,  rdummy, rdummy, rdummy, rdummy

    IF ( zpsoilout /= 0) THEN
      DO k=1,nzsoil
        WRITE (flunit) ((zpsoil(i,j,k),i=1,nx),j=1,ny)
      END DO
    END IF

    IF ( tsoilout /= 0 ) THEN
      DO is=0,nstyp
        DO k=1,nzsoil
          WRITE (flunit) ((tsoil(i,j,k,is),i=1,nx),j=1,ny)
        END DO
      END DO
    END IF

    IF ( qsoilout /= 0 ) THEN
      DO is=0,nstyp
        DO k=1,nzsoil
          WRITE (flunit) ((qsoil(i,j,k,is),i=1,nx),j=1,ny)
        END DO
      END DO
    END IF

    IF ( wcanpout /= 0 ) THEN
      DO is=0,nstyp
        WRITE (flunit) ((wetcanp(i,j,is),i=1,nx),j=1,ny)
      END DO
    END IF

    IF ( snowdout /= 0 ) THEN
      WRITE (flunit) ((snowdpth(i,j),i=1,nx),j=1,ny)
    END IF

    IF ( stypout /= 0 ) THEN
      DO is=1,nstyp
        WRITE (flunit) ((soiltyp(i,j,is),i=1,nx),j=1,ny)
      END DO
    ENDIF

    CLOSE ( flunit )
    CALL retunit ( flunit )

  ELSE IF (soildmp == 3) THEN

!-----------------------------------------------------------------------
!
!  Write out in HDF4.
!
!-----------------------------------------------------------------------

    intver = intver500  !  for the time being, in the future, we will
                        !  allow to dump data in the different version
                        !  intver will be assigned from input file
    IF (intver == intver500) THEN
      fmtver=fmtver500
    ELSE
      WRITE(6,'(/1x,a,i10,a/)')                                         &
          'Data format, intver=',intver,', not found. The Job stopped.'
      CALL arpsstop('arpstop called from wrtsoil.',1)
    END IF

    ALLOCATE(var3d (nxout, nyout, MAX(nzsoil,nstyp+1)), STAT = istat)
    ALLOCATE(var3di(nxout, nyout, nstyp),               STAT = istat)
    ALLOCATE(var4d (nxout, nyout, nzsoil,nstyp+1),      STAT = istat)

    DO jpy = 1, npyout
      DO ipx = 1, npxout

        ia = (ipx-1)*(nxout-3)
        ja = (jpy-1)*(nyout-3)

        IF (splitsoil > 0 .OR. mp_opt > 0 ) THEN
          CALL gtsplitfn(soiloutfl(1:lenbase),npxout,npyout,            &
                         loc_x,loc_y,ipx,jpy,                           &
                         0,0,0,lvldbg,outfilename,istat)
        ELSE
          outfilename = soiloutfl
          CALL fnversn(outfilename, lenbase)
        END IF

        CALL hdfopen(TRIM(outfilename), 2, sd_id)
        IF (sd_id < 0) THEN
          WRITE (6,*) "WRTSOIL: ERROR creating HDF4 file: ",            &
                      TRIM(outfilename)
          CALL arpsstop("ERROR creating soil HDF4 file.",1)
        END IF

        CALL hdfwrtc(sd_id,40,"fmtver",fmtver,  istat)
        CALL hdfwrti(sd_id, 'nstyp',   nstyp,   istat)
        CALL hdfwrti(sd_id, 'nx',      nxout,   istat)
        CALL hdfwrti(sd_id, 'ny',      nyout,   istat)
        CALL hdfwrti(sd_id, 'nzsoil',  nzsoil,  istat)
        CALL hdfwrtr(sd_id, 'dx',      dx,      istat)
        CALL hdfwrtr(sd_id, 'dy',      dy,      istat)
        CALL hdfwrti(sd_id, 'mapproj', mapproj, istat)
        CALL hdfwrtr(sd_id, 'trulat1', trulat1, istat)
        CALL hdfwrtr(sd_id, 'trulat2', trulat2, istat)
        CALL hdfwrtr(sd_id, 'trulon', trulon, istat)
        CALL hdfwrtr(sd_id, 'sclfct', sclfct, istat)
        CALL hdfwrtr(sd_id, 'ctrlat', ctrlat, istat)
        CALL hdfwrtr(sd_id, 'ctrlon', ctrlon, istat)

        nstyp = MAX(1,nstyp)

        IF ( zpsoilout /= 0) THEN

          DO k = 1,nzsoil
            DO j = 1, nyout
              DO i = 1, nxout
                var3d(i,j,k) = zpsoil(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL hdfwrt3d(var3d,nxout,nyout,nzsoil,sd_id,0,0,             &
                        'zpsoil','Soil level depth','m',                &
                        itmp,atmp1,atmp2)
        END IF

        IF ( tsoilout /= 0 ) THEN

          DO is = 1, nstyp+1
            DO k = 1,nzsoil
              DO j = 1, nyout
                DO i = 1, nxout
                  var4d(i,j,k,is) = tsoil(i+ia,j+ja,k,is-1)
                END DO
              END DO
            END DO
          END DO

          CALL hdfwrt4d(var4d,nxout,nyout,nzsoil,nstyp+1,sd_id,0,0,     &
                        'tsoil','Soil temperature','K',                 &
                        itmp,atmp1,atmp2)
        END IF

        IF ( qsoilout /= 0 ) THEN

          DO is = 1, nstyp+1
            DO k = 1,nzsoil
              DO j = 1, nyout
                DO i = 1, nxout
                  var4d(i,j,k,is) = qsoil(i+ia,j+ja,k,is-1)
                END DO
              END DO
            END DO
          END DO

          CALL hdfwrt4d(var4d,nxout,nyout,nzsoil,nstyp+1,sd_id,0,0,     &
                        'qsoil','Soil moisture','m3/m3',                &
                        itmp,atmp1,atmp2)
        END IF

        IF ( wcanpout /= 0 ) THEN

          DO is = 1, nstyp+1
            DO j = 1, nyout
              DO i = 1, nxout
                var3d(i,j,is) = wetcanp(i+ia,j+ja,is-1)
              END DO
            END DO
          END DO

          CALL hdfwrt3d(var3d,nxout,nyout,nstyp+1,sd_id,0,0,            &
                        'wetcanp','Canopy water amount','fraction',     &
                        itmp,atmp1,atmp2)
        END IF

        IF ( snowdout /= 0 ) THEN
          DO j = 1, nyout
            DO i = 1,nxout
              var3d(i,j,1) = snowdpth(i+ia,j+ja)
            END DO
          END DO
          CALL hdfwrt2d(var3d,nxout,nyout,sd_id,0,0,                    &
                        'snowdpth','Snow depth','m',itmp)
        END IF

        IF ( stypout /= 0 ) THEN
          DO k = 1,nstyp
            DO j = 1, nyout
              DO i = 1, nxout
                var3di(i,j,k) = soiltyp(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL hdfwrt3di(var3di,nxout,nyout,nstyp,sd_id,0,0,            &
                         'soiltyp','Soil type','index')
        ENDIF

        CALL hdfclose(sd_id,istat)
        IF (istat /= 0) THEN
          WRITE (6,*) "WRTSOIL: ERROR on closing file ",trim(soiloutfl), &
                      " (status",istat,")"
        END IF
      END DO
    END DO

    DEALLOCATE( var3d, var3di, var4d )

  ELSE IF (soildmp == 7) THEN
!-----------------------------------------------------------------------
!
!  Write out in NetCDF 3.0
!
!-----------------------------------------------------------------------

    nstyp=max(1,nstyp)

    ALLOCATE(var3d (nxout-1,nyout-1,MAX(nzsoil,nstyp+1)), STAT = istat)
    ALLOCATE(var3di(nxout-1,nyout-1,nstyp),               STAT = istat)
    ALLOCATE(var4d (nxout-1,nyout-1,nzsoil,nstyp+1),      STAT = istat)

!-----------------------------------------------------------------------
!
!  Define soil file dimension and variables
!
!-----------------------------------------------------------------------

    DO jpy = 1, npyout
      DO ipx = 1, npxout

        ia = (ipx-1)*(nxout-3)
        ja = (jpy-1)*(nyout-3)

        IF (splitsoil > 0 .OR. mp_opt > 0 ) THEN
          CALL gtsplitfn(soiloutfl(1:lenbase),npxout,npyout,            &
                         loc_x,loc_y,ipx,jpy,                           &
                         0,0,0,lvldbg,outfilename,istat)
        ELSE
          outfilename = soiloutfl
          CALL fnversn(outfilename, lenbase)
        END IF

        CALL netopen(TRIM(outfilename), 'C', sd_id)

        CALL net_define_soil(sd_id,nxout,nyout,nzsoil,nstyp, dx,dy,     &
                             mapproj,sclfct,                            &
                             trulat1,trulat2,trulon,ctrlat,ctrlon,      &
                             zpsoilout,tsoilout,qsoilout,wcanpout,      &
                             snowdout,stypout,istat)

        IF ( zpsoilout /= 0) THEN
          DO k = 1,nzsoil
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3d(i,j,k) = zpsoil(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3d(sd_id,0,0,'ZPSOIL',var3d,nxout-1,nyout-1,nzsoil)
        END IF

        IF ( tsoilout /= 0 ) THEN
          DO is = 0, nstyp
            DO k = 1,nzsoil
              DO j = 1,nyout-1
                DO i = 1,nxout-1
                  var4d(i,j,k,is+1) = tsoil(i+ia,j+ja,k,is)
                END DO
              END DO
            END DO
          END DO
          CALL netwrt4d(sd_id,0,0,'TSOIL',var4d,nxout-1,nyout-1,nzsoil,nstyp+1)
        END IF

        IF ( qsoilout /= 0 ) THEN
          DO is = 0, nstyp
            DO k = 1,nzsoil
              DO j = 1,nyout-1
                DO i = 1,nxout-1
                  var4d(i,j,k,is+1) = qsoil(i+ia,j+ja,k,is)
                END DO
              END DO
            END DO
          END DO
          CALL netwrt4d(sd_id,0,0,'QSOIL',var4d,nxout-1,nyout-1,nzsoil,nstyp+1)
        END IF

        IF ( wcanpout /= 0 ) THEN
          DO is = 0, nstyp
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3d(i,j,is+1) = wetcanp(i+ia,j+ja,is)
              END DO
            END DO
          END DO
          CALL netwrt3d(sd_id,0,0,'WETCANP',var3d,nxout-1,nyout-1,nstyp+1)
        END IF

        IF ( snowdout /= 0 ) THEN
          DO j = 1,nyout-1
            DO i = 1,nxout-1
              var3d(i,j,1) = snowdpth(i+ia,j+ja)
            END DO
          END DO
          CALL netwrt2d(sd_id,0,0,'SNOWDPTH',var3d,nxout-1,nyout-1)
        END IF

        IF ( stypout /= 0 ) THEN
          DO k = 1,nstyp
            DO j = 1,nyout-1
              DO i = 1,nxout-1
                var3di(i,j,k) = soiltyp(i+ia,j+ja,k)
              END DO
            END DO
          END DO
          CALL netwrt3di(sd_id,0,0,'SOILTYP',var3di,nxout-1,nyout-1,nstyp)
        ENDIF

        CALL netclose(sd_id)
      END DO
    END DO

    DEALLOCATE(var3d,var4d,var3di)

  ELSE

    ! alternate dump format ...
    WRITE(6,*) 'The supported soil data dump format are ',              &
               'binary (soildmp=1), HDF4 no compressed (soildmp = 3), ',&
               'and NetCDF (soildmp = 7).'
    CALL arpsstop('Soil data dump format is not supported.',1)

  END IF

  RETURN
END SUBROUTINE wrtsoil
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTJOINSOIL                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtjoinsoil( nx,ny,nzsoil,nstyps, soiloutfl, dx,dy, zpsoil,  &
                  mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,  &
                  zpsoilout,tsoilout,qsoilout,wcanpout,snowdout,        &
                  tsoil,qsoil,wetcanp,snowdpth,soiltyp )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Join the soil model variables and write into file soiloutfl.
!
!  NOTE:
!
!    Any changes made in this subroutine should also be made in wrtsoil above.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/20/2003
!  Modified from subroutine wrtsoil. Subroutine wrtsoil and wrfjoinsoil
!  should be changed simultaneously.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx        Number of model grid points in the x-dir. (east/west)
!  ny        Number of model grid points in the y-dir. (north/south)
!  nzsoil    Number of model grid points in the soil.
!
!  soiloutfl Name of the soil model data output file
!
!  zpsoilout Flag for output of soil level depth
!  tsoilout  Flag for output of soil temperature
!  qsoilout  Flag for output of soil moisture
!  wcanpout  Flag for output of Canopy water amount
!
!  tsoil     Soil temperature (K)
!  qsoil     Soil moisture (m3/m3)
!  wetcanp   Canopy water amount
!  snowdpth  Snow depth (m)
!  soiltyp   Soil type in model domain
!
!  OUTPUT:
!
!
!  Output written to the surface data file, soilinfl:
!
!  nx       Number of model grid points in the x-direction
!  ny       Number of model grid points in the y-direction
!  nzsoil   Number of model grid points in the soil
!
!  mapproj  Type of map projection used to setup the analysis grid.
!  trulat1  1st real true latitude of map projection.
!  trulat2  2nd real true latitude of map projection.
!  trulon   Real true longitude of map projection.
!  sclfct   Map scale factor. At latitude = trulat1 and trulat2
!
!  dx       Model grid spacing in the x-direction east-west (meters)
!  dy       Model grid spacing in the y-direction east-west (meters)
!  ctrlat   Lat. at the origin of the model grid (deg. N)
!  ctrlon   Lon. at the origin of the model grid (deg. E)
!
!  tsoilout Flag for output of soil temperature
!  qsoilout Flag for output of soil moisture
!  wcanpout Flag for output of Canopy water amount
!
!  tsoil    Soil temperature (K)
!  qsoil    Soil moisture (m3/m3)
!  wetcanp  Canopy water amount
!  snowdpth Snow depth (m)
!  soiltyp  Soil type in model domain
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx            ! Number of grid points in the x-direction
  INTEGER :: ny            ! Number of grid points in the y-direction
  INTEGER :: nzsoil        ! Number of grid points in the soil

  INTEGER :: nstyps        ! Number of soil types at each grid point
  CHARACTER(LEN=*) :: soiloutfl ! Name of the output file

  REAL :: dx
  REAL :: dy
  INTEGER :: mapproj       ! Map projection
  REAL :: trulat1          ! 1st real true latitude of map projection
  REAL :: trulat2          ! 2nd real true latitude of map projection
  REAL :: trulon           ! Real true longitude of map projection
  REAL :: sclfct           ! Map scale factor
  REAL :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon           ! Center longitude of the model domain (deg. E)

  INTEGER :: zpsoilout     ! Flag for output of zpsoil
  INTEGER :: tsoilout      ! Flag for output of tsoil
  INTEGER :: qsoilout      ! Flag for output of qsoil
  INTEGER :: wcanpout      ! Flag for output of wetcanp
  INTEGER :: stypout       ! Flag for output of soiltyp (set to 1 if any of
                           ! the above flags are set)
  INTEGER :: snowdout      ! Flag for output of snowdpth

  REAL :: zpsoil (nx,ny,nzsoil)          ! Soil depth (m)
  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m3/m3)
  REAL :: wetcanp(nx,ny,0:nstyps)    ! Canopy water amount

  REAL :: snowdpth(nx,ny)            ! Snow depth (m)

  INTEGER :: soiltyp (nx,ny,nstyps)  ! Soil type in model domain
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,is
  INTEGER :: flunit, lfn
  INTEGER :: idummy
  REAL    :: rdummy
  INTEGER :: ierr

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
  REAL :: atmp1(1),atmp2(1)
         ! unused arrays in hdf routines since NO COMPRESSION

  INTEGER :: istat, sd_id

!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver500
  INTEGER  :: intver,intver500

  PARAMETER (fmtver500='* 005.00 GrADS Soilvar Data',intver500=500)

!
! MP variables
!
  INTEGER :: nxlg,nylg

  REAL,    ALLOCATABLE :: out4d(:,:,:,:), out3d(:,:,:)
  INTEGER, ALLOCATABLE :: soiltyplg(:,:,:)
  REAL,    ALLOCATABLE :: out2d(:,:)

  REAL,    ALLOCATABLE :: var3d (:,:,:)
  INTEGER, ALLOCATABLE :: var3di(:,:,:)
  REAL,    ALLOCATABLE :: var4d (:,:,:,:)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

  idummy = 0
  rdummy = 0.0

  IF (soildmp == 0) RETURN

  IF ( tsoilout+qsoilout+wcanpout /= 0 ) THEN
    stypout = 1
  ELSE
    stypout = 0
  END IF

!
! To avoid nstyp > nstyps, in such a case tsoil, qsoil, wetcanp will be
! undefined for soil types diemension great than nstyps. -- WYH.
!
  nstyp = MIN(nstyp, nstyps)

  ALLOCATE(out4d(nxlg,nylg,nzsoil,0:nstyps),      STAT = istat)
  ALLOCATE(out3d(nxlg,nylg,MAX(nzsoil,nstyps+1)), STAT = istat)
  ALLOCATE(out2d(nxlg,nylg),                      STAT = istat)
  ALLOCATE(soiltyplg(nxlg,nylg,nstyps),           STAT = istat)

  lfn = LEN_TRIM(soiloutfl)
  CALL fnversn(soiloutfl, lfn)

  IF ( myproc == 0 )     &
    WRITE (6,'(1x,/,2a)') 'WRTJOINSOIL: Opening file ',TRIM(soiloutfl)

!-----------------------------------------------------------------------
!
!  Write out in Fortran unformatted.
!
!-----------------------------------------------------------------------

  IF (soildmp == 1) THEN

    IF(myproc == 0) THEN
      CALL getunit( flunit )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(soiloutfl, '-F f77 -N ieee', ierr)

      intver = intver500  !  for the time being, in the future, we will
                          !  allow to dump data in the different version
                          !  intver will be assigned from input file
      IF (intver == intver500) THEN
        fmtver=fmtver500
      ELSE
        WRITE(6,'(/1x,a,i10,a/)')                        &
          'Data format, intver=',intver,', not found. The Job stopped.'
        CALL arpsstop('arpstop called from wrtjoinsoil.',1)
      END IF

      OPEN (UNIT=flunit, FILE=trim(soiloutfl),STATUS='unknown',        &
            FORM='unformatted', ACCESS='sequential')

      WRITE (flunit) fmtver

      WRITE (flunit) nxlg, nylg, nzsoil

      WRITE (flunit) mapproj,tsoilout, qsoilout,                       &
                     wcanpout,      0,snowdout,  stypout, zpsoilout,   &
                     idummy,   idummy,  idummy,   idummy,  idummy,     &
                     idummy,   idummy,  idummy,   idummy,  nstyp

      WRITE (flunit)    dx,     dy, ctrlat, ctrlon,trulat1,            &
                   trulat2, trulon, sclfct, rdummy, rdummy,            &
                   rdummy,  rdummy, rdummy, rdummy, rdummy,            &
                   rdummy,  rdummy, rdummy, rdummy, rdummy
    END IF

    IF ( zpsoilout /= 0) THEN
      CALL mpimerge3d(zpsoil,nx,ny,nzsoil,out3d)
      IF(myproc == 0) THEN
        DO k=1,nzsoil
          WRITE (flunit) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
        END DO
      END IF
    END IF

    IF ( tsoilout /= 0 ) THEN
      CALL mpimerge4d(tsoil,nx,ny,nzsoil,nstyps+1,out4d)
      IF(myproc == 0) THEN
        DO is=0,nstyp
          DO k=1,nzsoil
            WRITE (flunit) ((out4d(i,j,k,is),i=1,nxlg),j=1,nylg)
          END DO
        END DO
      END IF
    END IF

    IF ( qsoilout /= 0 ) THEN
      CALL mpimerge4d(qsoil,nx,ny,nzsoil,nstyps+1,out4d)
      IF(myproc == 0) THEN
        DO is=0,nstyp
          DO k=1,nzsoil
            WRITE (flunit) ((out4d(i,j,k,is),i=1,nxlg),j=1,nylg)
          END DO
        END DO
      END IF
    END IF

    IF ( wcanpout /= 0 ) THEN
      CALL mpimerge3d(wetcanp, nx,ny,nstyps+1,out3d)
      IF(myproc == 0) THEN
        DO is=1,nstyp+1
          WRITE (flunit) ((out3d(i,j,is),i=1,nxlg),j=1,nylg)
        END DO
      END IF
    END IF

    IF ( snowdout /= 0 ) THEN
      CALL mpimerge2d(snowdpth,nx,ny,out2d)
      IF(myproc == 0 ) WRITE (flunit) ((out2d(i,j),i=1,nxlg),j=1,nylg)
    END IF

    IF ( stypout /= 0 ) THEN
      CALL mpimerge3di(soiltyp,nx,ny,nstyps,soiltyplg)
      IF(myproc == 0) THEN
        DO is=1,nstyp
          WRITE (flunit) ((soiltyplg(i,j,is),i=1,nxlg),j=1,nylg)
        END DO
      END IF
    ENDIF

    IF (myproc == 0) THEN
      CLOSE ( flunit )
      CALL retunit ( flunit )
    END IF

  ELSE IF (soildmp == 3) THEN

!-----------------------------------------------------------------------
!
!  Write out in HDF4.
!
!-----------------------------------------------------------------------

    intver = intver500  !  for the time being, in the future, we will
                        !  allow to dump data in the different version
                        !  intver will be assigned from input file
    IF (intver == intver500) THEN
      fmtver=fmtver500
    ELSE
      WRITE(6,'(/1x,a,i10,a/)')                        &
          'Data format, intver=',intver,', not found. The Job stopped.'
      CALL arpsstop('arpstop called from wrtjoinsoil.',1)
    END IF

    IF(myproc == 0) THEN
      CALL hdfopen(trim(soiloutfl), 2, sd_id)
      IF (sd_id < 0) THEN
        WRITE (6,*) "WRTJOINSOIL: ERROR creating HDF4 file: ",              &
                    trim(soiloutfl)
        CALL arpsstop("ERROR creating soil HDF4 file.",1)
      END IF

      CALL hdfwrtc(sd_id,40,"fmtver",fmtver,istat)
      CALL hdfwrti(sd_id, 'nstyp',  nstyp,  istat)
      CALL hdfwrti(sd_id, 'nx',     nxlg,   istat)
      CALL hdfwrti(sd_id, 'ny',     nylg,   istat)
      CALL hdfwrti(sd_id, 'nzsoil', nzsoil, istat)
      CALL hdfwrtr(sd_id, 'dx', dx, istat)
      CALL hdfwrtr(sd_id, 'dy', dy, istat)
      CALL hdfwrti(sd_id, 'mapproj', mapproj, istat)
      CALL hdfwrtr(sd_id, 'trulat1', trulat1, istat)
      CALL hdfwrtr(sd_id, 'trulat2', trulat2, istat)
      CALL hdfwrtr(sd_id, 'trulon',  trulon,  istat)
      CALL hdfwrtr(sd_id, 'sclfct',  sclfct,  istat)
      CALL hdfwrtr(sd_id, 'ctrlat',  ctrlat,  istat)
      CALL hdfwrtr(sd_id, 'ctrlon',  ctrlon,  istat)
    END IF

    nstyp=max(1,nstyp)

    IF ( zpsoilout /= 0) THEN
      CALL mpimerge3d(zpsoil,nx,ny,nzsoil,out3d)
      IF(myproc == 0)     &
        CALL hdfwrt3d(out3d,nxlg,nylg,nzsoil,sd_id,0,0,       &
                      'zpsoil','Soil level depth','m',        &
                      itmp,atmp1,atmp2)
    END IF

    IF ( tsoilout /= 0 ) THEN
      CALL mpimerge4d(tsoil,nx,ny,nzsoil,nstyps+1,out4d)
      IF(myproc == 0)    &
        CALL hdfwrt4d(out4d,nxlg,nylg,nzsoil,nstyp+1,sd_id,0,0, &
                      'tsoil','Soil temperature','K',         &
                      itmp,atmp1,atmp2)
    END IF

    IF ( qsoilout /= 0 ) THEN
      CALL mpimerge4d(qsoil,nx,ny,nzsoil,nstyps+1,out4d)
      IF(myproc == 0)    &
        CALL hdfwrt4d(out4d,nxlg,nylg,nzsoil,nstyp+1,sd_id,0,0, &
                      'qsoil','Soil moisture','m3/m3',        &
                      itmp,atmp1,atmp2)
    END IF

    IF ( wcanpout /= 0 ) THEN
      CALL mpimerge3d(wetcanp, nx,ny,nstyps+1,out3d)
      IF(myproc == 0)   &
        CALL hdfwrt3d(out3d,nxlg,nylg,nstyp+1,sd_id,0,0,          &
                    'wetcanp','Canopy water amount','fraction',   &
                    itmp,atmp1,atmp2)
    END IF

    IF ( snowdout /= 0 ) THEN
      CALL mpimerge2d(snowdpth,nx,ny,out2d)
      IF(myproc == 0) CALL hdfwrt2d(out2d,nxlg,nylg,sd_id,0,0,   &
                                 'snowdpth','Snow depth','m',itmp)
    END IF

    IF ( stypout /= 0 ) THEN
      CALL mpimerge3di(soiltyp,nx,ny,nstyps,soiltyplg)
      IF(myproc == 0) CALL hdfwrt3di(soiltyplg,nxlg,nylg,nstyp,  &
                                     sd_id,0,0,'soiltyp',        &
                                     'Soil type','index')
    ENDIF

    IF(myproc == 0) THEN
      CALL hdfclose(sd_id,istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "WRTJOINSOIL: ERROR on closing file ",trim(soiloutfl),    &
                    " (status",istat,")"
      END IF
    END IF

  ELSE IF (soildmp == 7) THEN
!-----------------------------------------------------------------------
!
!  Write out in NetCDF 3.0
!
!-----------------------------------------------------------------------

    nstyp=max(1,nstyp)

    ALLOCATE(var3d (nxlg-1,nylg-1,MAX(nzsoil,nstyp+1)), STAT = istat)
    ALLOCATE(var3di(nxlg-1,nylg-1,nstyp),               STAT = istat)
    ALLOCATE(var4d (nxlg-1,nylg-1,nzsoil,nstyp+1),      STAT = istat)

    IF (myproc == 0) THEN

!-----------------------------------------------------------------------
!
!  Define soil file dimension and variables
!
!-----------------------------------------------------------------------

      CALL netopen(TRIM(soiloutfl), 'C', sd_id)

      CALL net_define_soil(sd_id,nxlg,nylg,nzsoil,nstyp,dx,dy,mapproj,  &
           sclfct,trulat1,trulat2,trulon,ctrlat,ctrlon,                 &
           zpsoilout,tsoilout,qsoilout,wcanpout,snowdout,stypout,istat)

    END IF

    IF ( zpsoilout /= 0) THEN
      CALL mpimerge3d(zpsoil,nx,ny,nzsoil,out3d)
      IF(myproc == 0) THEN
        DO k = 1,nzsoil
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3d(sd_id,0,0,'ZPSOIL',var3d,nxlg-1,nylg-1,nzsoil)
      END IF
    END IF

    IF ( tsoilout /= 0 ) THEN
      CALL mpimerge4d(tsoil,nx,ny,nzsoil,nstyps+1,out4d)
      IF(myproc == 0) THEN
        DO is = 1,nstyp+1
          DO k = 1,nzsoil
            DO j = 1,nylg-1
              DO i = 1,nxlg-1
                var4d(i,j,k,is) = out4d(i,j,k,is-1)
              END DO
            END DO
          END DO
        END DO
        CALL netwrt4d(sd_id,0,0,'TSOIL',var4d,nxlg-1,nylg-1,nzsoil,nstyp+1)
      END IF
    END IF

    IF ( qsoilout /= 0 ) THEN
      CALL mpimerge4d(qsoil,nx,ny,nzsoil,nstyps+1,out4d)
      IF(myproc == 0) THEN
        DO is = 1,nstyp+1
          DO k = 1,nzsoil
            DO j = 1,nylg-1
              DO i = 1,nxlg-1
                var4d(i,j,k,is) = out4d(i,j,k,is-1)
              END DO
            END DO
          END DO
        END DO
       CALL netwrt4d(sd_id,0,0,'QSOIL',var4d,nxlg-1,nylg-1,nzsoil,nstyp+1)
      END IF
    END IF

    IF ( wcanpout /= 0 ) THEN
      CALL mpimerge3d(wetcanp,nx,ny,nstyps+1,out3d)
      IF(myproc == 0) THEN
        DO is = 1,nstyp+1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,is) = out3d(i,j,is)
            END DO
          END DO
        END DO
        CALL netwrt3d(sd_id,0,0,'WETCANP',var3d,nxlg-1,nylg-1,nstyp+1)
      END IF
    END IF

    IF ( snowdout /= 0 ) THEN
      CALL mpimerge2d(snowdpth,nx,ny,out2d)
      IF(myproc == 0) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = out2d(i,j)
          END DO
        END DO
        CALL netwrt2d(sd_id,0,0,'SNOWDPTH',var3d,nxlg-1,nylg-1)
      END IF
    END IF

    IF ( stypout /= 0 ) THEN
      CALL mpimerge3di(soiltyp,nx,ny,nstyp,soiltyplg)
      IF(myproc == 0) THEN
        DO k = 1,nstyp
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3di(i,j,k) = soiltyplg(i,j,k)
            END DO
          END DO
        END DO
        CALL netwrt3di(sd_id,0,0,'SOILTYP',var3di,nxlg-1,nylg-1,nstyp)
      END IF
    END IF

    IF(myproc == 0) CALL netclose(sd_id)

    DEALLOCATE(var3d,var3di)
    DEALLOCATE(var4d)

  ELSE

    ! alternate dump format ...
    WRITE(6,*) 'The supported soil data dump format are ',              &
               'binary (soildmp=1), HDF4 no compressed (soildmp = 3),', &
               'and NetCDF (soildmp = 7).'
    CALL arpsstop('Soil data dump format is not supported.',1)

  END IF

  DEALLOCATE(out4d, out3d, out2d, soiltyplg)

  RETURN
END SUBROUTINE wrtjoinsoil
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE READSOIL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE readsoil( nx,ny,nzsoil,nstyps,soilfile,dx,dy,zpsoil,         &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,            &
                   tsoil,qsoil,wetcanp,snowdpth,soiltyp )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the soil variables from file soilfile.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  04/20/95
!
!  MODIFICATION HISTORY:
!
!  2000/03/29 (Gene Bassett)
!  Removed the globcst.inc include.
!
!  2000/03/29 (Gene Bassett)
!  Added HDF4 format.
!
!  2000/10/27 (Gene Bassett)
!  Added soiltyp to soil file to track consistency of soil data.
!
!  05/14/2002 (J. Brotzge)
!  Added additional arrays, changed call statements to allow for
!  multiple soil schemes.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of model grid points in the x-dir. (east/west)
!  ny       Number of model grid points in the y-dir. (north/south)
!  nzsoil   Number of model grid points in the soil.
!
!  mapproj  Type of map projection used to setup the analysis grid.
!  trulat1  1st real true latitude of map projection.
!  trulat2  2nd real true latitude of map projection.
!  trulon   Real true longitude of map projection.
!  sclfct   Map scale factor. At latitude = trulat1 and trulat2
!
!  dx       Model grid spacing in the x-direction east-west (meters)
!  dy       Model grid spacing in the y-direction east-west (meters)
!  ctrlat    Lat. at the origin of the model grid (deg. N)
!  ctrlon    Lon. at the origin of the model grid (deg. E)
!
!  OUTPUT:
!
!  tsoil    Soil temperature (K)
!  qsoil    Soil moisture (m3/m3)
!  wetcanp  Canopy water amount
!  snowdpth Snow depth (m)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx            ! Number of grid points in the x-direction
  INTEGER :: ny            ! Number of grid points in the y-direction
  INTEGER :: nzsoil        ! Number of grid points in the soil.
  INTEGER :: nstyps        ! Number of soil types for each grid point

  CHARACTER (LEN=*) :: soilfile ! Name of the soil file

  REAL :: dx, dy
  INTEGER :: mapproj       ! Map projection
  REAL :: trulat1          ! 1st real true latitude of map projection
  REAL :: trulat2          ! 2nd real true latitude of map projection
  REAL :: trulon           ! Real true longitude of map projection
  REAL :: sclfct           ! Map scale factor
  REAL :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon           ! Center longitude of the model domain (deg. E)

  INTEGER :: zpsoilin, tsoilin, qsoilin, wcanpin
  INTEGER :: snowdin,  snowcin, stypin
  INTEGER :: tsfcin,   wsfcin,  wdpin       ! for backward compatibility Zuwen He, 07/01/02

  REAL    :: zpsoil (nx,ny,nzsoil)          ! Soil depths (m)
  REAL    :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL    :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m3/m3)
  REAL    :: wetcanp(nx,ny,0:nstyps)        ! Canopy water amount
  INTEGER :: soiltyp (nx,ny,nstyps)         ! Soil type in model domain

  REAL    :: snowdpth(nx,ny)                ! Snow depth (m)

  REAL,    ALLOCATABLE :: zpsoil_in (:,:,:)   ! Soil level depth (m)
  REAL,    ALLOCATABLE :: tsoil_in  (:,:,:,:) ! Soil temperature (K)
  REAL,    ALLOCATABLE :: qsoil_in  (:,:,:,:) ! Soil moisture (m3/m3)
  REAL,    ALLOCATABLE :: wetcanp_in(:,:,:)   ! Canopy water amount
  INTEGER, ALLOCATABLE :: soiltyp_in(:,:,:)   ! Soil type in model domain
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxin,nyin,nzsoilin, nstypin
  REAL    :: dxin,dyin
  INTEGER :: mprojin
  REAL    :: trlat1in, trlat2in, trlonin, sclfctin
  REAL    :: ctrlonin, ctrlatin

  INTEGER :: flunit
  INTEGER :: idummy, nstyp1
  REAL    :: rdummy

  INTEGER :: i,j,k,is
  INTEGER :: istat, ierr

  INTEGER :: ireturn

  CHARACTER (LEN=256) :: savename            !Message passing code.
  CHARACTER :: amm*2, ayear*4, aday*2

  INTEGER(KIND=selected_int_kind(4)) :: itmp(nx,ny,nzsoil,0:nstyps)
  REAL :: atmp1(nzsoil),atmp2(nzsoil)

  INTEGER :: sd_id
!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER(LEN=40) :: fmtver,fmtver410,fmtver500
  INTEGER           :: intver,intver410,intver500

  PARAMETER (fmtver410='* 004.10 GrADS Soilvar Data',intver410=410)
  PARAMETER (fmtver500='* 005.00 GrADS Soilvar Data',intver500=500)

  CHARACTER (LEN=40) :: fmtverin

  REAL, allocatable :: tem1(:,:,:) ! Temporary array

  REAL,    ALLOCATABLE :: var3d (:,:,:)
  INTEGER, ALLOCATABLE :: var3di(:,:,:)
  REAL,    ALLOCATABLE :: var4d (:,:,:,:)
!
!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
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
!  Open the surface data file. Read the parameters first, check if
!  the data are consistant with the model. If everything is OK, then
!  read the surface data, soiltyp, vegtyp, lai, and roufns.
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    savename(1:256) = soilfile(1:256)
    CALL gtsplitfn(savename,1,1,loc_x,loc_y,1,1,0,0,1,lvldbg,soilfile,ireturn)
  END IF

  IF (myproc == 0) WRITE (6,'(1x,/,3a)') "READSOIL: reading in ",       &
                 "external supplied soil data from file ",trim(soilfile)

!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------

  IF (soilfmt == 1) THEN

!-----------------------------------------------------------------------
!
!  Fortran unformatted dump.
!
!-----------------------------------------------------------------------

    CALL getunit( flunit )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(soilfile, '-F f77 -N ieee', ierr)

    OPEN (UNIT=flunit,FILE=trim(soilfile),FORM='unformatted',           &
        STATUS='old',IOSTAT=istat)

    IF( istat /= 0 ) THEN

      WRITE(6,'(/1x,a,i2,/1x,a/)')                                      &
          'Error occured when opening the surface data file '           &
          //soilfile//' using FORTRAN unit ',flunit,                    &
          ' Program stopped in READSOIL.'
      CALL arpsstop("arpsstop called from READSOIL while opening file",1)

    END IF

    WRITE(6,'(/1x,a,/1x,a,i2/)')                                        &
        'This run will start from an external supplied soil ',          &
        'data file '//soilfile//' using FORTRAN unit ',flunit

    READ (flunit,ERR=997) fmtverin
!
! 07/01/2002 Zuwen He
!
! The following code is not a safe practice.
!
! However, this may be the only way to distinguish versions
! prior to 500.
!
    IF (fmtverin == fmtver500) THEN
      intver=intver500
    ELSE
      WRITE(6,'(/1x,a/)')   &
          'WARNING: Incoming data format are older than version 5.00!!! '
    END IF

    WRITE(6,'(/1x,a,a/)') 'Incoming data format, fmtverin=',fmtverin

    READ (flunit,ERR=998) nxin,nyin,nzsoilin

    GOTO 996

    997 WRITE(6,'(/1x,a,a/)')                        &
      'Incoming data format: fmtver=fmtver410. Data read-in may be wrong.'

    CLOSE (flunit)
    OPEN (UNIT=flunit,FILE=trim(soilfile),FORM='unformatted',           &
        STATUS='old',IOSTAT=istat)

    READ (flunit,ERR=998) nxin,nyin
    nzsoilin=2
    intver=intver410   ! there is no fmtverin prior to version 500
    fmtver=fmtver410
    WRITE(6,'(/1x,a/,a/)')   &
          'WARNING: Incoming data format are to read as version 4.10'

    996 CONTINUE

    IF (intver == intver410) THEN

      READ (flunit,ERR=998) mprojin,tsfcin,tsoilin,wsfcin,wdpin,      &
                        wcanpin,snowcin,snowdin,stypin,zpsoilin,        &
                        idummy, idummy, idummy, idummy,idummy,          &
                        idummy, idummy, idummy, idummy,nstypin

    ELSE IF (intver >= intver500) THEN

      READ (flunit,ERR=998) mprojin,tsoilin,qsoilin,                    &
                        wcanpin,snowcin,snowdin,stypin,zpsoilin,        &
                        idummy, idummy, idummy, idummy,idummy,          &
                        idummy, idummy, idummy, idummy,nstypin

    END IF

    READ (flunit,ERR=998) dxin,dyin, ctrlatin,ctrlonin,trlat1in,        &
                         trlat2in,trlonin,sclfctin,rdummy,rdummy,       &
                         rdummy,rdummy,rdummy,rdummy,rdummy,            &
                         rdummy,rdummy,rdummy,rdummy,rdummy

  ELSE IF (soilfmt == 3) THEN     !HDF4

!-----------------------------------------------------------------------
!
!  HDF4 format.
!
!-----------------------------------------------------------------------

    CALL hdfopen(trim(soilfile), 1, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "READSOIL: ERROR opening ",                           &
                 trim(soilfile)," for reading."
      GO TO 998
    END IF

    CALL hdfrdc(sd_id,40,"fmtver",fmtverin,istat)

!
! 07/01/2002  Zuwen He
!
! The following code is a dangerous practice
! but it may be the only way to distinguish
! versions prior to 500.
!
    IF (fmtverin == fmtver500) THEN
      intver=intver500
    ELSE
      intver=intver410  ! prior to 500, there is no fmtver variable
      istat=0
      WRITE(6,'(/1x,a/,a/)')   &
          'WARNING: Incoming data format are older than version 5.00!!! ', &
          'It is to be read as if it version 4.10!!! '
    END IF

    CALL hdfrdi(sd_id,"nstyp",nstypin,istat)
    CALL hdfrdi(sd_id,"nx",   nxin,   istat)
    CALL hdfrdi(sd_id,"ny",   nyin,   istat)

    IF (intver >= intver500) THEN
      CALL hdfrdi(sd_id,"nzsoil",nzsoilin,istat)
    ELSE
      nzsoilin = 2  ! prior to version 500, it is 2 layer soil
    END IF

    CALL hdfrdr(sd_id,"dx",dxin,istat)
    CALL hdfrdr(sd_id,"dy",dyin,istat)
    CALL hdfrdi(sd_id,"mapproj",mprojin, istat)
    CALL hdfrdr(sd_id,"trulat1",trlat1in,istat)
    CALL hdfrdr(sd_id,"trulat2",trlat2in,istat)
    CALL hdfrdr(sd_id,"trulon", trlonin, istat)
    CALL hdfrdr(sd_id,"sclfct", sclfctin,istat)
    CALL hdfrdr(sd_id,"ctrlat", ctrlatin,istat)
    CALL hdfrdr(sd_id,"ctrlon", ctrlonin,istat)

  ELSE IF (soilfmt == 7) THEN     ! NetCDF

!-----------------------------------------------------------------------
!
!  NetCDF format.
!
!-----------------------------------------------------------------------

    CALL netopen(TRIM(soilfile), 'R', sd_id)

    CALL net_get_soil(sd_id,nxin,nyin,nzsoilin,nstypin,dxin,dyin,       &
          mprojin,sclfctin,trlat1in,trlat2in,trlonin,ctrlatin,ctrlonin, &
          zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,stypin,istat)

    intver = intver500

  ELSE IF (sfcfmt == 2) THEN  !Read data directly
    nstypin = 1

  ELSE

    ! alternate dump format ...
    WRITE(6,'(1x,3a)') 'The supported soil data format are ',           &
              'binary (soilfmt=1), HDF4 no compressed (soilfmt = 3). ', &
              'and NetCDF (soilfmt = 7).'
    CALL arpsstop('Soil data format is not supported.',1)

  END IF

  nstyp1 = MAX( nstypin, 1 )

  ALLOCATE (zpsoil_in(nx,ny,nzsoil),stat=istat)
  CALL check_alloc_status(istat, "readsoil:zpsoil_in")

  ALLOCATE (tsoil_in(nx,ny,nzsoil,0:nstyp1),stat=istat)
  CALL check_alloc_status(istat, "readsoil:tsoil_in")

  ALLOCATE (qsoil_in(nx,ny,nzsoil,0:nstyp1),stat=istat)
  CALL check_alloc_status(istat, "readsoil:qsoil_in")

  ALLOCATE (wetcanp_in(nx,ny,0:nstyp1),stat=istat)
  CALL check_alloc_status(istat, "readsoil:wetcanp_in")

  ALLOCATE (soiltyp_in(nx,ny,nstyp1),stat=istat)
  CALL check_alloc_status(istat, "readsoil:soiltyp_in")

!-----------------------------------------------------------------------
!
!  Check the data file for consistent grid parameters.
!
!-----------------------------------------------------------------------

  IF (soilfmt /= 2) THEN                !(OASIS testing, sfcfmt = 2)

    CALL checkgrid2d(nx,ny,nxin,nyin,                                   &
              dx,dy,ctrlat,ctrlon,                                      &
              mapproj,trulat1,trulat2,trulon,sclfct,                    &
              dxin,dyin,ctrlatin,ctrlonin,                              &
              mprojin,trlat1in,trlat2in,trlonin,sclfctin,ireturn)

    IF (ireturn /= 0) THEN
      WRITE (6,*) "ERROR: READSOIL, grid parameter mismatch"
      CALL arpsstop("arpsstop called from READSOIL, parameter mismatch",1)
    END IF

    IF (myproc == 0) WRITE (6,'(a,a//a,i2/6(a,e15.8/))')                &
      ' The map projection and griding information for the ',           &
      ' surface data: ',                                                &
              ' Projection:                 ', mprojin,                 &
              ' The 1st real true latitude: ', trlat1in,                &
              ' The 2nd real true latitude: ', trlat2in,                &
              ' The real true longitude:    ', trlonin,                 &
              ' Map scale factor:           ', sclfctin,                &
              ' Latitude  at the origin:    ', ctrlatin,                &
              ' Longitude at the origin:    ', ctrlonin

    IF (nzsoil /= nzsoilin) THEN
      WRITE(6,*)                                                        &
            'ERROR: nzsoilin in soil file not equal to expected nzsoil.'
      WRITE(6,*) 'nzsoil = ',nzsoil,' nzsoilin = ',nzsoilin
      CALL arpsstop("arpsstop called from READSOIL, parameter mismatch",1)
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in the soil data from the soil data file.
!
!-----------------------------------------------------------------------
!
  IF (intver == intver410) THEN
    ALLOCATE (tem1(nx,ny,0:nstyp1),stat=istat)  ! for reading old version
                                                ! tsfc,tsoil,wetsfc,wetdp
  END IF

  IF (soilfmt == 1) THEN    ! Fortran unformatted

   IF (intver == intver410) THEN

       WRITE (6, '(a/8(a,i2/))')  &
                    ' Surface data flags for: ',                        &
                    '        zpsoilin --', zpsoilin,                    &
                    '        tsfcin   --', tsfcin,                      &
                    '        tsoilin  --', tsoilin,                     &
                    '        wsfcin   --', wsfcin,                      &
                    '        wdpin    --', wdpin,                       &
                    '        wcanpin  --', wcanpin,                     &
                    '        snowdin  --', snowdin,                     &
                    '        stypin   --', stypin

     ELSE IF (intver == intver500) THEN

       WRITE (6, '(a/6(a,i2/))')  &
                    ' Surface data flags for: ',                        &
                    '        zpsoilin --', zpsoilin,                    &
                    '        tsoilin  --', tsoilin,                     &
                    '        qsoilin  --', qsoilin,                     &
                    '        wcanpin  --', wcanpin,                     &
                    '        snowdin  --', snowdin,                     &
                    '        stypin   --', stypin

    END IF

    IF (intver == intver410) THEN

      zpsoil_in(:,:,1)=0.
      zpsoil_in(:,:,2)=-1.

    ELSE IF (intver >= intver500) THEN

      IF (zpsoilin /= 0) THEN
        WRITE(6,'(a)') 'Read in the soil depth data'
        DO k=1,nzsoilin
          READ (flunit,ERR=998) ((zpsoil_in(i,j,k),i=1,nxin),j=1,nyin)
        END DO
      END IF

    END IF  ! intver

    IF (intver == intver410) THEN

      IF ( tsfcin /= 0 ) THEN
        WRITE(6,'(a)') 'Read in the surface skin temperature data'
        IF ( nstyp1 == 1 ) THEN
          READ (flunit,ERR=998) ((tem1(i,j,0),i=1,nxin),j=1,nyin)
          tsoil_in(:,:,1,0)=tem1(:,:,0)
        ELSE
          READ (flunit,ERR=998) tem1
          tsoil_in(:,:,1,:)=tem1(:,:,:)
        END IF
      ELSE
        WRITE(6,'(a)') 'Variable tsfc is not in the data set.'
      END IF

      IF ( tsoilin /= 0 ) THEN
        WRITE(6,'(a)') 'Read in the deep soil temperature data'
        IF ( nstyp1 == 1 ) THEN
          READ (flunit,ERR=998) ((tem1(i,j,0),i=1,nxin),j=1,nyin)
          tsoil_in(:,:,2,0)=tem1(:,:,0)
        ELSE
          READ (flunit,ERR=998) tem1
          tsoil_in(:,:,2,:)=tem1(:,:,:)
        END IF
      ELSE
        WRITE(6,'(a)') 'Variable tsoil is not in the data set.'
      END IF

      IF ( wsfcin /= 0 ) THEN
        WRITE(6,'(a)') 'Read in the skin soil moisture data'
        IF ( nstyp1 == 1 ) THEN
          READ (flunit,ERR=998) ((tem1(i,j,0),i=1,nxin),j=1,nyin)
          qsoil_in(:,:,1,0)=tem1(:,:,0)
        ELSE
          READ (flunit,ERR=998) tem1
          qsoil_in(:,:,1,:)=tem1(:,:,:)
        END IF
      ELSE
      WRITE(6,'(a)') 'Variable wetsfc is not in the data set.'
      END IF

      IF ( wdpin /= 0 ) THEN
        WRITE(6,'(a)') 'Read in the deep soil moisture data'
        IF ( nstyp1 == 1 ) THEN
          READ (flunit,ERR=998) ((tem1(i,j,0),i=1,nxin),j=1,nyin)
          qsoil_in(:,:,2,0)=tem1(:,:,0)
        ELSE
          READ (flunit,ERR=998) tem1
          qsoil_in(:,:,2,:)=tem1(:,:,:)
        END IF
      ELSE
        WRITE(6,'(a)') 'Variable wetdp is not in the data set.'
      END IF

      IF ( wcanpin /= 0 ) THEN
        IF ( nstyp1 == 1 ) THEN
          READ (flunit,ERR=998) ((wetcanp_in(i,j,0),i=1,nxin),j=1,nyin)
        ELSE
          READ (flunit,ERR=998) wetcanp_in
        END IF
      ELSE
        WRITE (6, '(a)') 'Variable wetcanp is not in the data set.'
      END IF

      IF ( snowcin /= 0 ) THEN
        WRITE (6, '(a)') 'File contains snowcvr -- discarding'
        READ (flunit,ERR=998)
      END IF

      IF ( snowdin /= 0 ) THEN
        WRITE (6, '(a)') 'Read in the snow depth data'
        READ (flunit,ERR=998) snowdpth
      ELSE
        WRITE (6, '(a)') 'Variable snowdpth is not in the data set.'
      END IF

      IF ( stypin /= 0 ) THEN
        WRITE (6, '(a)') 'Read soil type of soil data.'
        READ (flunit,ERR=998) soiltyp_in
      END IF

    ELSE IF (intver >= intver500) THEN

      IF ( tsoilin /= 0 ) THEN
        DO is=0,nstypin
          WRITE(6,'(a,i4)') 'Read in the soil temperature for soil type ',is
          DO k=1,nzsoilin
            READ (flunit,ERR=998) ((tsoil_in(i,j,k,is),i=1,nxin),j=1,nyin)
          END DO
        END DO
      ELSE
        WRITE(6,'(a)') 'Variable tsoil is not in the data set.'
      END IF

      IF ( qsoilin /= 0 ) THEN
        DO is=0,nstypin
          WRITE(6,'(a,i4)') 'Read in the soil moisture data for soil type ',is
          DO k=1,nzsoilin
            READ (flunit,ERR=998) ((qsoil_in(i,j,k,is),i=1,nxin),j=1,nyin)
          END DO
        END DO
      ELSE
        WRITE(6,'(a)') 'Variable qsoil is not in the data set.'
      END IF

      IF ( wcanpin /= 0 ) THEN
        DO is=0,nstypin
          WRITE (6, '(a,i4)') 'Read in the canopy water amount data for soil type ',is
          READ (flunit,ERR=998) ((wetcanp_in(i,j,is),i=1,nxin),j=1,nyin)
        END DO
      ELSE
        WRITE (6, '(a)') 'Variable wetcanp is not in the data set.'
      END IF

      IF ( snowcin /= 0 ) THEN
        WRITE (6, '(a)') 'File contains snowcvr -- discarding'
        READ (flunit,ERR=998)
      END IF

      IF ( snowdin /= 0 ) THEN
        WRITE (6, '(a)') 'Read in the snow depth data'
        READ (flunit,ERR=998) ((snowdpth(i,j),i=1,nxin),j=1,nyin)
      ELSE
        WRITE (6, '(a)') 'Variable snowdpth is not in the data set.'
      END IF

      IF ( stypin /= 0 ) THEN
        DO is=1,nstypin
          WRITE (6, '(a,i4)') 'Read soil type of soil data for soil type ',is
          READ (flunit,ERR=998) ((soiltyp_in(i,j,is),i=1,nxin),j=1,nyin)
        END DO
      END IF

    END IF

    CLOSE ( flunit )
    CALL retunit ( flunit )

  ELSE IF (soilfmt == 3) THEN

    IF (intver <= intver410) THEN

      WRITE(6,'( 1x,a)') 'WARNING: No zpsoil is defined in this version. '
      WRITE(6,'(10x,a)') 'Assume zpsoil_in(,,1)=0 and zpsoil(,,2)=-1.'
      zpsoil_in(:,:,1)=0.
      zpsoil_in(:,:,2)=-1.

      CALL hdfrd3d(sd_id,"tsfc",nxin,nyin,nstyp1+1,tem1,istat,          &
                   itmp,atmp1,atmp2)
      IF (istat > 1) GO TO 998
      IF (istat == 0) THEN
        WRITE(6,'(1x,a)') 'Read in the surface skin temperature data'
        tsfcin = 1
      ELSE
        WRITE(6,'(1x,a)') 'Variable tsfc is not in the data set.'
        tsfcin = 0
      END IF
      tsoil_in(:,:,1,:)=tem1(:,:,:)

      CALL hdfrd3d(sd_id,"tsoil",nxin,nyin,nstyp1+1,tem1,istat,         &
                   itmp,atmp1,atmp2)
      IF (istat > 1) GO TO 998
      IF (istat == 0) THEN
        WRITE(6,'(1x,a)') 'Read in the deep soil temperature data'
        tsoilin = 1
      ELSE
        WRITE(6,'(1x,a)') 'Variable tsoil is not in the data set.'
        tsoilin = 0
      END IF
      tsoil_in(:,:,2,:)=tem1(:,:,:)

      CALL hdfrd3d(sd_id,"wetsfc",nxin,nyin,nstyp1+1,tem1,istat,        &
                   itmp,atmp1,atmp2)
      IF (istat > 1) GO TO 998
      IF (istat == 0) THEN
        WRITE(6,'(1x,a)') 'Read in the skin soil moisture data'
        wsfcin = 1
      ELSE
        WRITE(6,'(1x,a)') 'Variable wetsfc is not in the data set.'
        wsfcin = 0
      END IF
      qsoil_in(:,:,1,:)=tem1(:,:,:)

      CALL hdfrd3d(sd_id,"wetdp",nxin,nyin,nstyp1+1,tem1,istat,         &
                   itmp,atmp1,atmp2)
      IF (istat > 1) GO TO 998
      IF (istat == 0) THEN
        WRITE(6,'(1x,a)') 'Read in the deep soil moisture data'
        wdpin = 1
      ELSE
        WRITE(6,'(1x,a)') 'Variable wetdp is not in the data set.'
        wdpin = 0
      END IF
      qsoil_in(:,:,2,:)=tem1(:,:,:)

    ELSE IF (intver >= intver500) THEN

      CALL hdfrd3d(sd_id,"zpsoil",nxin,nyin,nzsoil,zpsoil_in,istat,     &
                 itmp,atmp1,atmp2)
      IF (istat > 1) GO TO 998
      IF (istat == 0) THEN
        IF (myproc == 0) WRITE(6,'(1x,a)') 'Read in the soil layer depth data'
        zpsoilin = 1
      ELSE
        WRITE(6,'(1x,a)') 'Variable zpsoil is not in the data set.'
        zpsoilin = 0
      END IF

     CALL hdfrd4d(sd_id,"tsoil",nxin,nyin,nzsoil,nstyp1+1,tsoil_in,istat,&
                 itmp,atmp1,atmp2)
      IF (istat > 1) GO TO 998
      IF (istat == 0) THEN
        IF (myproc == 0) WRITE(6,'(1x,a)') 'Read in the soil temperature data'
        tsoilin = 1
      ELSE
        WRITE(6,'(1x,a)') 'Variable tsoil is not in the data set.'
        tsoilin = 0
      END IF

     CALL hdfrd4d(sd_id,"qsoil",nxin,nyin,nzsoil,nstyp1+1,qsoil_in,istat,   &
                   itmp,atmp1,atmp2)
      IF (istat > 1) GO TO 998
      IF (istat == 0) THEN
        IF (myproc == 0) WRITE(6,'(1x,a)') 'Read in the soil moisture data'
        qsoilin = 1
      ELSE
        WRITE(6,'(1x,a)') 'Variable qsoil is not in the data set.'
        qsoilin = 0
      END IF

    END IF

    CALL hdfrd3d(sd_id,"wetcanp",nxin,nyin,nstyp1+1,wetcanp_in,istat,   &
                 itmp,atmp1,atmp2)
    IF (istat > 1) GO TO 998
    IF (istat == 0) THEN
      IF (myproc == 0) WRITE (6, '(1x,a)') 'Read in the canopy water amount data'
      wcanpin = 1
    ELSE
      WRITE (6, '(1x,a)') 'Variable wetcanp is not in the data set.'
      wcanpin = 0
    END IF

    CALL hdfrd2d(sd_id,"snowdpth",nxin,nyin,snowdpth,istat,itmp)
    IF (istat > 1) GO TO 998
    IF (istat == 0) THEN
      IF (myproc == 0) WRITE (6, '(1x,a)') 'Read in the snow depth data'
      snowdin = 1
    ELSE
      WRITE (6, '(1x,a)') 'Variable snowdpth is not in the data set.'
      snowdin = 0
    END IF

    CALL hdfrd3di(sd_id,"soiltyp",nxin,nyin,nstyp1,soiltyp_in,istat)
    IF (istat > 1) GO TO 998
    IF (istat == 0) THEN
      IF (myproc == 0) WRITE (6, '(1x,a)') 'Read soil type of soil data'
      stypin = 1
    ELSE
      WRITE (6, '(1x,a)') 'Soil type of soil data is not in the data set.'
      stypin = 0
    END IF

    CALL hdfclose(sd_id, istat)

  ELSE IF (soilfmt == 7) THEN

    ALLOCATE(var3d (nxin-1,nyin-1,MAX(nzsoil,nstyp1+1)), STAT = istat)
    ALLOCATE(var3di(nxin-1,nyin-1,nstyp1),               STAT = istat)
    ALLOCATE(var4d (nxin-1,nyin-1,nzsoil,nstyp1+1),      STAT = istat)

    IF (zpsoilin == 1) THEN
      CALL netread3d(sd_id,0,0,"ZPSOIL",nxin-1,nyin-1,nzsoil,var3d)
      DO k = 1, nzsoil
        DO j = 1,nyin-1
          DO i = 1,nxin-1
            zpsoil_in(i,j,k) = var3d(i,j,k)
          END DO
        END DO
      END DO
      CALL edgfill(zpsoil_in,1,nxin,1,nxin-1,1,nyin,1,nyin-1,1,nzsoil,1,nzsoil)
    END IF

    IF (tsoilin == 1) THEN
      CALL netread4d(sd_id,0,0,"TSOIL",nxin-1,nyin-1,nzsoil,nstyp1+1,var4d)
      DO is = 0, nstyp1
        DO k = 1, nzsoil
          DO j = 1,nyin-1
            DO i = 1,nxin-1
              tsoil_in(i,j,k,is) = var4d(i,j,k,is+1)
            END DO
          END DO
        END DO
        CALL edgfill(tsoil_in(:,:,:,is),1,nxin,1,nxin-1,1,nyin,1,nyin-1,1,nzsoil,1,nzsoil)
      END DO
    END IF

    IF (qsoilin == 1) THEN
      CALL netread4d(sd_id,0,0,"QSOIL",nxin-1,nyin-1,nzsoil,nstyp1+1,var4d)
      DO is = 0, nstyp1
        DO k = 1, nzsoil
          DO j = 1,nyin-1
            DO i = 1,nxin-1
              qsoil_in(i,j,k,is) = var4d(i,j,k,is+1)
            END DO
          END DO
        END DO
        CALL edgfill(qsoil_in(:,:,:,is),1,nxin,1,nxin-1,1,nyin,1,nyin-1,1,nzsoil,1,nzsoil)
      END DO
    END IF

    IF (wcanpin == 1) THEN
      CALL netread3d(sd_id,0,0,"WETCANP",nxin-1,nyin-1,nstyp1+1,var3d)
      DO is = 0, nstyp1
        DO j = 1,nyin-1
          DO i = 1,nxin-1
            wetcanp_in(i,j,is) = var3d(i,j,is+1)
          END DO
        END DO
      END DO
      CALL edgfill(wetcanp_in,1,nxin,1,nxin-1,1,nyin,1,nyin-1,1,nstyp1+1,1,nstyp1+1)
    END IF

    IF (snowdin == 1) THEN
      CALL netread2d(sd_id,0,0,"SNOWDPTH",nxin-1,nyin-1,var3d)
      DO j = 1,nyin-1
        DO i = 1,nxin-1
          snowdpth(i,j) = var3d(i,j,1)
        END DO
      END DO
      CALL edgfill(snowdpth,1,nxin,1,nxin-1,1,nyin,1,nyin-1,1,1,1,1)
    END IF

    IF (stypin == 1) THEN
      CALL netread3di(sd_id,0,0,"SOILTYP",nxin-1,nyin-1,nstyp1,var3di)
      DO is = 1, nstyp1
        DO j = 1,nyin-1
          DO i = 1,nxin-1
            soiltyp_in(i,j,is) = var3di(i,j,is)
          END DO
        END DO
      END DO
      CALL iedgfill(soiltyp_in,1,nxin,1,nxin-1,1,nyin,1,nyin-1,1,nstyp1,1,nstyp1)
    END IF

    CALL netclose(sd_id)

    DEALLOCATE(var3d,var3di)
    DEALLOCATE(var4d)

    ! alternate dump format ...

  ELSE IF (soilfmt == 2) THEN     ! Data read directly (JAB)
                                  ! OASIS code
    CALL initztime(ayear,amm,aday)

    CALL readjsoil(nx,ny,nzsoil, nstyp,ayear,amm,aday,ztime,              &
        zpsoil,tsoil,qsoil,wetcanp,snowdpth)

  END IF

!-----------------------------------------------------------------------
!
! Consistency adjustment
!
!-----------------------------------------------------------------------

  IF (stypin == 0) THEN

    WRITE (6,'(1x,2a)')  "READSOIL: WARNING, no check made for ",       &
         "consistency between soil types in surface and soil data sets."

    nstypin = MIN(nstyp1, nstyps)

    IF (zpsoilin == 1) zpsoil (:,:,:) = zpsoil_in(:,:,:)
    IF (tsoilin == 1)  tsoil  (:,:,:,0:nstypin) = tsoil_in  (:,:,:,0:nstypin)
    IF (qsoilin == 1)  qsoil  (:,:,:,0:nstypin) = qsoil_in  (:,:,:,0:nstypin)
    IF (wcanpin == 1)  wetcanp(:,:,0:nstypin)   = wetcanp_in(:,:,0:nstypin)

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp1,nstyp,tsoil,qsoil,wetcanp)

  ELSE

    IF (nstyp1 == 1) THEN
      tsoil_in  (:,:,:,1) = tsoil_in  (:,:,:,0)
      qsoil_in  (:,:,:,1) = qsoil_in  (:,:,:,0)
      wetcanp_in(:,:,1)   = wetcanp_in(:,:,0)
    END IF

    IF (soilmodel_forced == 0) THEN

      CALL remap_soil_vars(nx,ny,nzsoil,nstyp1,nstyp,                   &
                tsoil_in,qsoil_in,wetcanp_in,soiltyp_in,                &
                tsfcin,tsoilin,wsfcin,wdpin,qsoilin,wcanpin,            &
                intver,                                                 &
                tsoil,qsoil,wetcanp,soiltyp)

    END IF

  END IF

  IF (mp_opt > 0) soilfile(1:256) = savename(1:256)

  IF (intver == intver410) DEALLOCATE (tem1)
  ! for reading old version tsfc,tsoil,wetsfc,wetdp

  DEALLOCATE (tsoil_in,stat=istat)
  DEALLOCATE (qsoil_in,stat=istat)
  DEALLOCATE (wetcanp_in,stat=istat)
  DEALLOCATE (soiltyp_in,stat=istat)

  IF (intver == intver410) THEN

    IF (tsfcin /= tsoilin .OR. wsfcin /= wdpin) THEN

    WRITE (6,'(a,a/,a/,a/)')   &
          'READSOIL: WARNING: The soilvar data is of version ', fmtver410,  &
          '. The inconsistency flag between tsfcin and tsoilin, ',  &
          ' or between wsfin and wdpin, may cause some problems. '
    END IF

    tsoilin = max(tsfcin,tsoilin)
    qsoilin = max(wsfcin,wdpin)

  END IF

! Correct only for flint, otherwise flint will never end
!  RETURN
  GO TO 999

  998   WRITE (6,'(/a,i2/a)')                                           &
         '     Read error in soil data file '                           &
         //soilfile//' with the I/O unit ',flunit,                      &
         'The model will STOP in subroutine READSOIL.'

  CALL arpsstop("arpsstop called from READSOIL reading surface data",1)

  999 CONTINUE

  RETURN

END SUBROUTINE readsoil
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE READSPLITSOIL             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE readsplitsoil( nx,ny,nzsoil,nstyps,soilfile,dx,dy,zpsoil,    &
                mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,    &
                zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,               &
                tsoil,qsoil,wetcanp,snowdpth,soiltyp )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the soil variables from file soilfile. Split and scatter to
!  MP processors from the root process.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/30/2002
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx        Number of model grid points in the x-dir. (east/west)
!  ny        Number of model grid points in the y-dir. (north/south)
!  nzsoil    Number of model grid points in the soil.
!
!  mapproj   Type of map projection used to setup the analysis grid.
!  trulat1   1st real true latitude of map projection.
!  trulat2   2nd real true latitude of map projection.
!  trulon    Real true longitude of map projection.
!  sclfct    Map scale factor. At latitude = trulat1 and trulat2
!
!  dx        Model grid spacing in the x-direction east-west (meters)
!  dy        Model grid spacing in the y-direction east-west (meters)
!  ctrlat    Lat. at the origin of the model grid (deg. N)
!  ctrlon    Lon. at the origin of the model grid (deg. E)
!
!  OUTPUT:
!
!  tsoil     Soil temperature (K)
!  qsoil     Soil moisture (m3/m3)
!  wetcanp   Canopy water amount
!  snowdpth  Snow depth (m)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx            ! Number of grid points in the x-direction
  INTEGER :: ny            ! Number of grid points in the y-direction
  INTEGER :: nzsoil        ! Number of grid points in the soil.
  INTEGER :: nstyps        ! Number of soil types for each grid point

  CHARACTER (LEN=*) :: soilfile ! Name of the soil file

  REAL :: dx
  REAL :: dy
  INTEGER :: mapproj       ! Map projection
  REAL :: trulat1          ! 1st real true latitude of map projection
  REAL :: trulat2          ! 2nd real true latitude of map projection
  REAL :: trulon           ! Real true longitude of map projection
  REAL :: sclfct           ! Map scale factor
  REAL :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon           ! Center longitude of the model domain (deg. E)

  INTEGER :: wcanpin,  snowdin, snowcin, stypin
  INTEGER :: zpsoilin, tsoilin, qsoilin
  INTEGER :: tsfcin,   wsfcin,  wdpin
                           ! for backward compatibility Zuwen He, 07/01/02

  REAL :: zpsoil (nx,ny,nzsoil)          ! Soil depths (m)
  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m3/m3)
  REAL :: wetcanp(nx,ny,0:nstyps)        ! Canopy water amount
  INTEGER :: soiltyp (nx,ny,nstyps)      ! Soil type in model domain

  REAL :: snowdpth(nx,ny)                ! Snow depth (m)

  REAL,    ALLOCATABLE :: zpsoil_in (:,:,:)     ! Soil level depth (m)
  REAL,    ALLOCATABLE :: tsoil_in  (:,:,:,:)   ! Soil temperature (K)
  REAL,    ALLOCATABLE :: qsoil_in  (:,:,:,:)   ! Soil moisture (m3/m3)
  REAL,    ALLOCATABLE :: wetcanp_in(:,:,:)     ! Canopy water amount
  INTEGER, ALLOCATABLE :: soiltyp_in(:,:,:)     ! Soil type in model domain
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxin,nyin,nzsoilin,nstypin
  INTEGER :: mprojin
  REAL    :: dxin,dyin
  REAL    :: trlat1in, trlat2in, trlonin, sclfctin
  REAL    :: ctrlonin, ctrlatin

  INTEGER :: flunit
  INTEGER :: idummy, nstyp1
  REAL    :: rdummy

  INTEGER :: i,j,k,is
  INTEGER :: istat, ierr

  INTEGER :: ireturn

  CHARACTER :: amm*2, ayear*4, aday*2

  INTEGER(KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:,:)
  REAL :: atmp1(nzsoil),atmp2(nzsoil)
         ! unused arrays in hdf routines since NO COMPRESSION

  INTEGER :: sd_id

  INTEGER              :: nxlg, nylg, n3rd
  INTEGER, ALLOCATABLE :: var3di(:,:,:)
  REAL,    AlLOCATABLE :: var3d(:,:,:), var4d(:,:,:,:)

  INTEGER, ALLOCATABLE :: net3di(:,:,:)
  REAL,    AlLOCATABLE :: net3d(:,:,:), net4d(:,:,:,:)
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver410,fmtver500
  INTEGER  :: intver,intver410,intver500

  PARAMETER (fmtver410='* 004.10 GrADS Soilvar Data',intver410=410)
  PARAMETER (fmtver500='* 005.00 GrADS Soilvar Data',intver500=500)

  CHARACTER (LEN=40) :: fmtverin

  REAL, ALLOCATABLE  :: tem1(:,:,:) ! Temporary array

!
!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3
!
!-----------------------------------------------------------------------
!
!  Open the surface data file. Read the parameters first, check if
!  the data are consistant with the model. If everything is OK, then
!  read the surface data, soiltyp, vegtyp, lai, and roufns.
!
!-----------------------------------------------------------------------
!

  IF (myproc == 0) THEN

    WRITE (6,'(1x,/,3a)') 'READSPLITSOIL: reading in external supplied',&
                          ' soil data from file - ',trim(soilfile)

!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------

    IF (soilfmt == 1) THEN

!-----------------------------------------------------------------------
!
!  Fortran unformatted dump.
!
!-----------------------------------------------------------------------

      CALL getunit( flunit )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(soilfile, '-F f77 -N ieee', ierr)

      OPEN(UNIT=flunit,FILE=trim(soilfile),FORM='unformatted',          &
           STATUS='old',IOSTAT=istat)

      IF( istat /= 0 ) THEN

        WRITE(6,'(/1x,a,i2,/1x,a/)')                                    &
            'Error occured when opening the surface data file '         &
            //soilfile//' using FORTRAN unit ',flunit,                  &
            ' Program stopped in READSPLITSOIL.'
        CALL arpsstop("arpsstop called from READSPLITSOIL while opening file",1)

      END IF

      WRITE(6,'(/1x,a,/1x,a,i2/)')                                      &
          'This run will start from an external supplied soil ',        &
          'data file '//soilfile//' using FORTRAN unit ',flunit

      READ (flunit,ERR=997) fmtverin
      !
      ! The following code is not a safe practice.
      !
      ! However, this may be the only way to distinguish versions
      ! prior to 500.
      !
      IF (fmtverin == fmtver500) THEN
        intver=intver500
      ELSE
        WRITE(6,'(/1x,a/)')   &
            'WARNING: Incoming data format are older than version 5.00!!! '
      END IF

      WRITE(6,'(/1x,a,a/)') 'Incoming data format, fmtverin=',fmtverin

      READ (flunit,ERR=998) nxin,nyin,nzsoilin

      GOTO 996

      997 WRITE(6,'(/1x,a,a/)')                        &
         'Incoming data format: fmtver=fmtver410. Data read-in may be wrong.'

      CLOSE (flunit)
      OPEN (UNIT=flunit,FILE=trim(soilfile),FORM='unformatted',           &
            STATUS='old',IOSTAT=istat)

      READ (flunit,ERR=998) nxin,nyin
      nzsoilin=2
      intver=intver410   ! there is no fmtverin prior to version 500
      fmtver=fmtver410
      WRITE(6,'(/1x,a/,a/)')   &
            'WARNING: Incoming data format are to read as version 4.10'

      996 CONTINUE

      IF (intver == intver410) THEN

        READ (flunit,ERR=998) mprojin,tsfcin,tsoilin,wsfcin,wdpin,      &
                          wcanpin,snowcin,snowdin,stypin,zpsoilin,        &
                          idummy, idummy, idummy, idummy,idummy,          &
                          idummy, idummy, idummy, idummy,nstypin

      ELSE IF (intver >= intver500) THEN

        READ (flunit,ERR=998) mprojin,tsoilin,qsoilin,                    &
                          wcanpin,snowcin,snowdin,stypin,zpsoilin,        &
                          idummy, idummy, idummy, idummy,idummy,          &
                          idummy, idummy, idummy, idummy,nstypin

      END IF

      READ (flunit,ERR=998) dxin,dyin, ctrlatin,ctrlonin,trlat1in,        &
                           trlat2in,trlonin,sclfctin,rdummy,rdummy,       &
                           rdummy,rdummy,rdummy,rdummy,rdummy,            &
                           rdummy,rdummy,rdummy,rdummy,rdummy

    ELSE IF (soilfmt == 3) THEN     !HDF4 format

      CALL hdfopen(trim(soilfile), 1, sd_id)
      IF (sd_id < 0) THEN
        WRITE (6,*) "READSPLITSOIL: ERROR opening ",                      &
                   trim(soilfile)," for reading."
        GO TO 998
      END IF

      CALL hdfrdc(sd_id,40,"fmtver",fmtverin,istat)
  !
  ! The following code is a dangerous practice
  ! but it may be the only way to distinguish
  ! versions prior to 500.
  !
      IF (fmtverin == fmtver500) THEN
        intver=intver500
      ELSE
        intver=intver410  ! prior to 500, there is no fmtver variable
        istat=0
        WRITE(6,'(/1x,a/,a/)')   &
            'WARNING: Incoming data format are older than version 5.00!!! ', &
            'It is to be read as if it version 4.10!!! '
      END IF

      CALL hdfrdi(sd_id,"nstyp",nstypin,istat)
      CALL hdfrdi(sd_id,"nx",   nxin,   istat)
      CALL hdfrdi(sd_id,"ny",   nyin,   istat)

      IF (intver >= intver500) THEN
        CALL hdfrdi(sd_id,"nzsoil",nzsoilin,istat)
      ELSE
        nzsoilin = 2  ! prior to version 500, it is 2 layer soil
      END IF

      CALL hdfrdr(sd_id,"dx",dxin,istat)
      CALL hdfrdr(sd_id,"dy",dyin,istat)
      CALL hdfrdi(sd_id,"mapproj",mprojin,istat)
      CALL hdfrdr(sd_id,"trulat1",trlat1in,istat)
      CALL hdfrdr(sd_id,"trulat2",trlat2in,istat)
      CALL hdfrdr(sd_id,"trulon",trlonin,istat)
      CALL hdfrdr(sd_id,"sclfct",sclfctin,istat)
      CALL hdfrdr(sd_id,"ctrlat",ctrlatin,istat)
      CALL hdfrdr(sd_id,"ctrlon",ctrlonin,istat)

    ELSE IF (soilfmt == 7) THEN     ! NetCDF 3.0 format

      CALL netopen(TRIM(soilfile), 'R', sd_id)

      CALL net_get_soil(sd_id,nxin,nyin,nzsoilin,nstypin,dxin,dyin,     &
          mprojin,sclfctin,trlat1in,trlat2in,trlonin,ctrlatin,ctrlonin, &
          zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,stypin,istat)

      intver = intver500

    ELSE IF (sfcfmt == 2) THEN  !Read data directly
      nstypin = 1

    ELSE

      ! alternate dump format ...
      WRITE(6,*) 'The supported soil data format are ',              &
               'binary (soilfmt=1) and HDF4 no compressed (soilfmt = 3).'
      CALL arpsstop('Soil data format is not supported.',1)

    END IF

    nxin = (nxin-3)/nproc_x + 3
    nyin = (nyin-3)/nproc_y + 3
    nstyp1 = MAX( nstypin, 1 )

  END IF   ! myproc == 0

  CALL mpupdatei(nstyp1, 1)
  CALL mpupdatei(intver, 1)
  CALL mpupdatei(nxin, 1)
  CALL mpupdatei(nyin, 1)
  CALL mpupdatei(nzsoilin, 1)
  CALL mpupdatei(nstypin, 1)
  CALL mpupdater(dxin, 1)
  CALL mpupdater(dyin, 1)
  CALL mpupdatei(mprojin, 1)
  CALL mpupdater(ctrlatin, 1)
  CALL mpupdater(ctrlonin, 1)
  CALL mpupdater(trlat1in, 1)
  CALL mpupdater(trlat2in, 1)
  CALL mpupdater(trlonin,  1)
  CALL mpupdater(sclfctin, 1)

  IF (soilfmt == 1 .OR. soilfmt == 7) THEN
    CALL mpupdatei(tsoilin, 1)
    CALL mpupdatei(wcanpin, 1)
    CALL mpupdatei(snowcin, 1)
    CALL mpupdatei(snowdin, 1)
    CALL mpupdatei(stypin,  1)
    CALL mpupdatei(zpsoilin,1)
    IF (intver == intver410) THEN
      CALL mpupdatei(tsfcin, 1)
      CALL mpupdatei(wsfcin, 1)
      CALL mpupdatei(wdpin, 1)
    ELSE
      CALL mpupdatei(qsoilin, 1)
    END IF
  END IF

  ALLOCATE (zpsoil_in(nx,ny,nzsoil),stat=istat)
  CALL check_alloc_status(istat,"readjoinsoil:zpsoil_in")

  ALLOCATE (tsoil_in(nx,ny,nzsoil,0:nstyp1),stat=istat)
  CALL check_alloc_status(istat,"readjoinsoil:tsoil_in")

  ALLOCATE (qsoil_in(nx,ny,nzsoil,0:nstyp1),stat=istat)
  CALL check_alloc_status(istat,"readjoinsoil:qsoil_in")

  ALLOCATE (wetcanp_in(nx,ny,0:nstyp1),stat=istat)
  CALL check_alloc_status(istat,"readjoinsoil:wetcanp_in")

  ALLOCATE (soiltyp_in(nx,ny,nstyp1),stat=istat)
  CALL check_alloc_status(istat,"readjoinsoil:soiltyp_in")

  IF (intver == intver410) THEN
    ALLOCATE (tem1(nx,ny,0:nstyp1),stat=istat)  ! for reading old version
                                                ! tsfc,tsoil,wetsfc,wetdp
  END IF

  ALLOCATE (itmp(nxlg, nylg, nzsoil, 0:nstyp1), stat = istat)
  CALL check_alloc_status(istat,"readjoinsoil:itmp")

  n3rd = MAX(nzsoil, nstyp1+1)

  ALLOCATE (var3d(nxlg, nylg, n3rd), stat = istat)
  CALL check_alloc_status(istat,"readjoinsoil:var3d")

  ALLOCATE (var3di(nxlg, nylg, nstyp1), stat = istat)
  CALL check_alloc_status(istat,"readjoinsoil:var3di")

  ALLOCATE (var4d(nxlg, nylg, nzsoil, nstyp1+1), stat = istat)
  CALL check_alloc_status(istat,"readjoinsoil:var4d")

!-----------------------------------------------------------------------
!
!  Check the data file for consistent grid parameters.
!
!-----------------------------------------------------------------------

  IF (soilfmt /= 2) THEN          !(OASIS testing, sfcfmt = 2)

    CALL checkgrid2d(nx,ny,nxin,nyin,                                   &
                   dx,dy,ctrlat,ctrlon,                                 &
                   mapproj,trulat1,trulat2,trulon,sclfct,               &
                   dxin,dyin,ctrlatin,ctrlonin,                         &
                   mprojin,trlat1in,trlat2in,trlonin,sclfctin,ireturn)

    IF (ireturn /= 0) THEN
      WRITE (6,*) "READSPLITSOIL: ERROR, grid parameter mismatch"
      CALL arpsstop("arpsstop called from READSPLITSOIL parameter mismatch",1)
    END IF

    IF (nzsoil /= nzsoilin) THEN
      WRITE(6,*)                                                          &
        'ERROR -- nzsoilin in soil file not equal to expected nzsoil.'
      WRITE(6,*)'nzsoil = ',nzsoil,' nzsoilin = ',nzsoilin
      CALL arpsstop("arpsstop called from READSPLITSOIL parameter mismatch",1)
    END IF

    IF (myproc == 0) WRITE (6,'(a,a//a,i2/6(a,e15.8/))')                &
               ' The map projection and griding information for the ',  &
               ' surface data: ',                                       &
               ' Projection:                 ', mprojin,                &
               ' The 1st real true latitude: ', trlat1in,               &
               ' The 2nd real true latitude: ', trlat2in,               &
               ' The real true longitude:    ', trlonin,                &
               ' Map scale factor:           ', sclfctin,               &
               ' Latitude  at the origin:    ', ctrlatin,               &
               ' Longitude at the origin:    ', ctrlonin

  END IF

!
!-----------------------------------------------------------------------
!
!  Read in the soil data from the soil data file.
!
!-----------------------------------------------------------------------
!
  IF (soilfmt == 1) THEN    ! Fortran unformatted

    IF (myproc == 0) THEN

      IF (intver == intver410) THEN

        WRITE (6, '(a/8(a,i2/))')      &
                        ' Surface data flags for: ',                    &
                        '        zpsoilin --', zpsoilin,                &
                        '        tsfcin   --', tsfcin,                  &
                        '        tsoilin  --', tsoilin,                 &
                        '        wsfcin   --', wsfcin,                  &
                        '        wdpin    --', wdpin,                   &
                        '        wcanpin  --', wcanpin,                 &
                        '        snowdin  --', snowdin,                 &
                        '        stypin   --', stypin
        var3d(:,:,1)=0.
        var3d(:,:,2)=-1.

      ELSE IF (intver == intver500) THEN

        WRITE (6, '(a/6(a,i2/))')      &
                        ' Surface data flags for: ',                    &
                        '        zpsoilin --', zpsoilin,                &
                        '        tsoilin  --', tsoilin,                 &
                        '        qsoilin  --', qsoilin,                 &
                        '        wcanpin  --', wcanpin,                 &
                        '        snowdin  --', snowdin,                 &
                        '        stypin   --', stypin

        IF (zpsoilin /= 0) THEN
          WRITE(6,'(a)') 'Read in the soil depth data'
          DO k=1,nzsoilin
            READ (flunit,ERR=998) ((var3d(i,j,k),i=1,nxlg),j=1,nylg)
          END DO
        END IF
     END IF  ! intver

   END IF  ! myproc == 0

   CALL mpisplit3d(var3d,nxin,nyin,nzsoilin,zpsoil_in)

   IF (intver == intver410) THEN

     IF (nstyp1 == 1) THEN
       n3rd = 1
     ELSE
       n3rd = nstyp1 + 1
     END IF

     IF ( tsfcin /= 0 ) THEN
       IF (myproc == 0) THEN
         WRITE(6,'(a)') 'Read in the surface skin temperature data'
         READ (flunit,ERR=998) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,n3rd)
       END IF
       CALL mpisplit3d(var3d, nxin,nyin,n3rd,tem1(:,:,0:n3rd-1))
       tsoil_in(:,:,1,0:n3rd-1)=tem1(:,:,0:n3rd-1)
     ELSE
       IF (myproc == 0) &
         WRITE(6,'(a)') 'Variable tsfc is not in the data set.'
     END IF

     IF ( tsoilin /= 0 ) THEN
       IF (myproc == 0) THEN
         WRITE(6,'(a)') 'Read in the deep soil temperature data'
         READ (flunit,ERR=998) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,n3rd)
       END IF
       CALL mpisplit3d(var3d, nxin,nyin,n3rd,tem1(:,:,0:n3rd-1))
       tsoil_in(:,:,2,0:n3rd-1)=tem1(:,:,0:n3rd-1)
     ELSE
       IF (myproc == 0)   &
         WRITE(6,'(a)') 'Variable tsoil is not in the data set.'
     END IF

     IF ( wsfcin /= 0 ) THEN
       IF (myproc == 0) THEN
         WRITE(6,'(a)') 'Read in the skin soil moisture data'
         READ (flunit,ERR=998) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,n3rd)
       END IF
       CALL mpisplit3d(var3d, nxin,nyin,n3rd,tem1(:,:,0:n3rd-1))
       qsoil_in(:,:,1,0:n3rd-1)=tem1(:,:,0:n3rd-1)
     ELSE
       IF (myproc == 0)   &
         WRITE(6,'(a)') 'Variable wetsfc is not in the data set.'
     END IF

     IF ( wdpin /= 0 ) THEN
       IF (myproc == 0) THEN
         WRITE(6,'(a)') 'Read in the deep soil moisture data'
         READ (flunit,ERR=998) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,n3rd)
       END IF
       CALL mpisplit3d(var3d, nxin,nyin,n3rd,tem1(:,:,0:n3rd-1))
       qsoil_in(:,:,2,0:n3rd-1)=tem1(:,:,0:n3rd-1)
     ELSE
       IF (myproc == 0)   &
         WRITE(6,'(a)') 'Variable wetdp is not in the data set.'
     END IF

     IF ( wcanpin /= 0 ) THEN
       IF (myproc == 0) THEN
         READ (flunit,ERR=998) (((var3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,n3rd)
       END IF
       CALL mpisplit3d(var3d, nxin,nyin,n3rd,wetcanp_in(:,:,0:n3rd-1))
     ELSE
       IF (myproc == 0)   &
         WRITE (6, '(a)') 'Variable wetcanp is not in the data set.'
     END IF

     IF ( snowcin /= 0 ) THEN
       IF (myproc == 0) THEN
         WRITE (6, '(a)') 'File contains snowcvr -- discarding'
         READ (flunit,ERR=998)
       END IF
     END IF

     IF ( snowdin /= 0 ) THEN
       IF (myproc == 0) THEN
         WRITE (6, '(a)') 'Read in the snow depth data'
         READ (flunit,ERR=998) ((var3d(i,j,1),i=1,nxlg),j=1,nylg)
       END IF
       CALL mpisplit3d(var3d, nxin,nyin,1,snowdpth)
     ELSE
       IF (myproc == 0)   &
         WRITE (6, '(a)') 'Variable snowdpth is not in the data set.'
     END IF

     IF ( stypin /= 0 ) THEN
       IF (myproc == 0) THEN
         WRITE (6, '(a)') 'Read soil type of soil data.'
         READ (flunit,ERR=998) (((var3di(i,j,k),i=1,nxlg),j=1,nylg),k=1,nstyp1)
       END IF
       CALL mpisplit3di(var3di, nxin,nyin,nstyp1,soiltyp_in)
     END IF

   ELSE IF (intver >= intver500) THEN

      IF ( tsoilin /= 0 ) THEN
        IF (myproc == 0) THEN
          DO is=0,nstyp1
            WRITE(6,'(a,i4)') 'Read in the soil temperature for soil type ',is
            DO k=1,nzsoilin
              READ (flunit,ERR=998) ((var4d(i,j,k,is+1),i=1,nxlg),j=1,nylg)
            END DO
          END DO
        END IF
        CALL mpisplit4d(var4d, nxin,nyin,nzsoil,nstyp1+1,tsoil_in)
      ELSE
        IF (myproc == 0)   &
          WRITE(6,'(a)') 'Variable tsoil is not in the data set.'
      END IF

      IF ( qsoilin /= 0 ) THEN
        IF (myproc == 0) THEN
          DO is=0,nstyp1
            WRITE(6,'(a,i4)') 'Read in the soil moisture data for soil type ',is
            DO k=1,nzsoilin
              READ (flunit,ERR=998) ((var4d(i,j,k,is+1),i=1,nxlg),j=1,nylg)
            END DO
          END DO
        END IF
        CALL mpisplit4d(var4d, nxin,nyin,nzsoil,nstyp1+1,qsoil_in)
      ELSE
        IF (myproc == 0)   &
          WRITE(6,'(a)') 'Variable qsoil is not in the data set.'
      END IF

      IF ( wcanpin /= 0 ) THEN
        IF (myproc == 0) THEN
          DO is=0,nstyp1
            WRITE (6, '(a,i4)') 'Read in the canopy water amount data for soil type ',is
            READ (flunit,ERR=998) ((var3d(i,j,is+1),i=1,nxlg),j=1,nylg)
          END DO
        END IF
        CALL mpisplit3d(var3d, nxin,nyin,nstyp1+1,wetcanp_in)
      ELSE
        IF (myproc == 0)   &
          WRITE (6, '(a)') 'Variable wetcanp is not in the data set.'
      END IF

      IF ( snowcin /= 0 ) THEN
        IF (myproc == 0) THEN
          WRITE (6, '(a)') 'File contains snowcvr -- discarding'
          READ (flunit,ERR=998)
        END IF
      END IF

      IF ( snowdin /= 0 ) THEN
        IF (myproc == 0) THEN
          WRITE (6, '(a)') 'Read in the snow depth data'
          READ (flunit,ERR=998) ((var3d(i,j,1),i=1,nxlg),j=1,nylg)
        END IF
        CALL mpisplit3d(var3d, nxin,nyin,1,snowdpth)
      ELSE
        IF (myproc == 0)   &
          WRITE (6, '(a)') 'Variable snowdpth is not in the data set.'
      END IF

      IF ( stypin /= 0 ) THEN
        IF (myproc == 0) THEN
          DO is=1,nstyp1
            WRITE (6, '(a,i4)') 'Read soil type of soil data for soil type ',is
            READ (flunit,ERR=998) ((var3d(i,j,is),i=1,nxlg),j=1,nylg)
          END DO
        END IF
        CALL mpisplit3di(var3di, nxin,nyin,nstyp1,soiltyp_in)
      END IF

    END IF

    IF (myproc == 0) THEN
      CLOSE ( flunit )
      CALL retunit ( flunit )
    END IF

  ELSE IF (soilfmt == 3) THEN

    IF (intver <= intver410) THEN

      IF (myproc == 0) THEN
        WRITE(6,'(a)') 'WARNING: No zpsoil is defined in this version. '
        WRITE(6,'(a)') 'Assume zpsoil_in(,,1)=0 and zpsoil(,,2)=-1.'
      END IF
      zpsoil_in(:,:,1)=0.
      zpsoil_in(:,:,2)=-1.

      IF (myproc == 0) THEN
        CALL hdfrd3d(sd_id,"tsfc",nxlg,nylg,nstyp1+1,var3d,istat,     &
                   itmp,atmp1,atmp2)
        IF (istat > 1) GO TO 998
        IF (istat == 0) THEN
          WRITE(6,'(a)') 'Read in the surface skin temperature data'
          tsfcin = 1
        ELSE
          WRITE(6,'(a)') 'Variable tsfc is not in the data set.'
          tsfcin = 0
        END IF
      END IF
      CALL mpupdatei(tsfcin, 1)
      CALL mpisplit3d(var3d, nxin, nyin,nstyp1+1,tem1)
      tsoil_in(:,:,1,:)=tem1(:,:,:)

      IF (myproc == 0) THEN
        CALL hdfrd3d(sd_id,"tsoil",nxlg,nylg,nstyp1+1,var3d,istat,    &
                   itmp,atmp1,atmp2)
        IF (istat > 1) GO TO 998
        IF (istat == 0) THEN
          WRITE(6,'(a)') 'Read in the deep soil temperature data'
          tsoilin = 1
        ELSE
          WRITE(6,'(a)') 'Variable tsoil is not in the data set.'
          tsoilin = 0
        END IF
      END IF
      CALL mpupdatei(tsoilin, 1)
      CALL mpisplit3d(var3d, nxin, nyin,nstyp1+1,tem1)
      tsoil_in(:,:,2,:)=tem1(:,:,:)

      IF (myproc == 0) THEN
        CALL hdfrd3d(sd_id,"wetsfc",nxlg,nylg,nstyp1+1,var3d,istat,   &
                   itmp,atmp1,atmp2)
        IF (istat > 1) GO TO 998
        IF (istat == 0) THEN
          WRITE(6,'(a)') 'Read in the skin soil moisture data'
          wsfcin = 1
        ELSE
          WRITE(6,'(a)') 'Variable wetsfc is not in the data set.'
          wsfcin = 0
        END IF
      END IF
      CALL mpupdatei(wsfcin, 1)
      CALL mpisplit3d(var3d, nxin, nyin,nstyp1+1,tem1)
      qsoil_in(:,:,1,:)=tem1(:,:,:)

      IF (myproc == 0) THEN
        CALL hdfrd3d(sd_id,"wetdp",nxlg,nylg,nstyp1+1,var3d,istat,                 &
                   itmp,atmp1,atmp2)
        IF (istat > 1) GO TO 998
        IF (istat == 0) THEN
          WRITE(6,'(a)') 'Read in the deep soil moisture data'
          wdpin = 1
        ELSE
          WRITE(6,'(a)') 'Variable wetdp is not in the data set.'
          wdpin = 0
        END IF
      END IF
      CALL mpupdatei(wdpin, 1)
      CALL mpisplit3d(var3d, nxin, nyin,nstyp1+1,tem1)
      qsoil_in(:,:,2,:)=tem1(:,:,:)

    ELSE IF (intver >= intver500) THEN

      IF (myproc == 0) THEN
        CALL hdfrd3d(sd_id,"zpsoil",nxlg,nylg,nzsoil,var3d,istat,         &
                 itmp,atmp1,atmp2)
        IF (istat > 1) GO TO 998
        IF (istat == 0) THEN
          WRITE(6,'(a)') 'Read in the soil layer depth data'
          zpsoilin = 1
        ELSE
          WRITE(6,'(a)') 'Variable zpsoil is not in the data set.'
          zpsoilin = 0
        END IF
      END IF
      CALL mpupdatei(zpsoilin, 1)
      CALL mpisplit3d(var3d, nxin, nyin,nzsoil,zpsoil_in)

      IF (myproc == 0) THEN
        CALL hdfrd4d(sd_id,"tsoil",nxlg,nylg,nzsoil,nstyp1+1,var4d,istat,   &
                 itmp,atmp1,atmp2)
        IF (istat > 1) GO TO 998
        IF (istat == 0) THEN
          WRITE(6,'(a)') 'Read in the soil temperature data'
          tsoilin = 1
        ELSE
          WRITE(6,'(a)') 'Variable tsoil is not in the data set.'
          tsoilin = 0
        END IF
      END IF
      CALL mpupdatei(tsoilin, 1)
      CALL mpisplit4d(var4d, nxin, nyin,nzsoil,nstyp1+1,tsoil_in)

      IF (myproc == 0) THEN
        CALL hdfrd4d(sd_id,"qsoil",nxlg,nylg,nzsoil,nstyp1+1,var4d,istat,   &
                   itmp,atmp1,atmp2)
        IF (istat > 1) GO TO 998
        IF (istat == 0) THEN
          WRITE(6,'(a)') 'Read in the soil moisture data'
          qsoilin = 1
        ELSE
          WRITE(6,'(a)') 'Variable qsoil is not in the data set.'
          qsoilin = 0
        END IF
      END IF
      CALL mpupdatei(qsoilin, 1)
      CALL mpisplit4d(var4d, nxin, nyin,nzsoil,nstyp1+1,qsoil_in)

    END IF

    IF (myproc == 0) THEN
      CALL hdfrd3d(sd_id,"wetcanp",nxlg,nylg,nstyp1+1,var3d,istat,         &
                 itmp,atmp1,atmp2)
      IF (istat > 1) GO TO 998
      IF (istat == 0) THEN
        WRITE (6, '(a)') 'Read in the canopy water amount data'
        wcanpin = 1
      ELSE
        WRITE (6, '(a)') 'Variable wetcanp is not in the data set.'
        wcanpin = 0
      END IF
    END IF
    CALL mpupdatei(wcanpin, 1)
    CALL mpisplit3d(var3d, nxin,nyin,nstyp1+1,wetcanp_in)

    IF (myproc == 0) THEN
      CALL hdfrd2d(sd_id,"snowdpth",nxlg,nylg,var3d,istat,itmp)
      IF (istat > 1) GO TO 998
      IF (istat == 0) THEN
        WRITE (6, '(a)') 'Read in the snow depth data'
        snowdin = 1
      ELSE
        WRITE (6, '(a)') 'Variable snowdpth is not in the data set.'
        snowdin = 0
      END IF
    END IF
    CALL mpupdatei(snowdin, 1)
    CALL mpisplit3d(var3d, nxin,nyin,1,snowdpth)

    IF (myproc == 0) THEN
      CALL hdfrd3di(sd_id,"soiltyp",nxlg,nylg,nstyp1,var3di,istat)
      IF (istat > 1) GO TO 998
      IF (istat == 0) THEN
        WRITE (6, '(a)') 'Read soil type of soil data'
        stypin = 1
      ELSE
        WRITE (6, '(a)') 'Soil type of soil data is not in the data set.'
        stypin = 0
      END IF
    END IF
    CALL mpupdatei(stypin, 1)
    CALL mpisplit3di(var3di, nxin,nyin,nstyp1,soiltyp_in)

    IF(myproc == 0) CALL hdfclose(sd_id, istat)

  ELSE IF (soilfmt == 7) THEN

    ALLOCATE(net3d(nxlg-1,nylg-1,MAX(nzsoil,nstyp1+1)), STAT = istat)
    ALLOCATE(net4d(nxlg-1,nylg-1,nzsoil,nstyp1+1),      STAT = istat)
    ALLOCATE(net3di(nxlg-1,nylg-1,nstyp1),              STAT = istat)

    IF (zpsoilin == 1) THEN
      IF (myproc == 0) THEN
        CALL netread3d(sd_id,0,0,"ZPSOIL",nxlg-1,nylg-1,nzsoil,net3d)
        WRITE(6,'(1x,a)') 'Read in the soil layer depth data.'
        DO k = 1,nzsoil
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,k) = net3d(i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nzsoil,1,nzsoil)
      END IF
      CALL mpisplit3d(var3d,nxin,nyin,nzsoil,zpsoil_in)
    END IF

    IF (tsoilin == 1) THEN
      IF (myproc == 0) THEN
        CALL netread4d(sd_id,0,0,"TSOIL",nxlg-1,nylg-1,nzsoil,nstyp1+1,net4d)
        WRITE(6,'(1x,a)') 'Read in the soil temperature data.'

        DO is = 1, nstyp1+1
          DO k = 1,nzsoil
            DO j = 1,nylg-1
              DO i = 1,nxlg-1
                var4d(i,j,k,is) = net4d(i,j,k,is)
              END DO
            END DO
          END DO
          CALL edgfill(var4d(:,:,:,is),1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nzsoil,1,nzsoil)
        END DO
      END IF
      CALL mpisplit4d(var4d,nxin,nyin,nzsoil,nstyp1+1,tsoil_in)
    END IF

    IF (qsoilin == 1) THEN
      IF (myproc == 0) THEN
        CALL netread4d(sd_id,0,0,"QSOIL",nxlg-1,nylg-1,nzsoil,nstyp1+1,net4d)
        WRITE(6,'(1x,a)') 'Read in the soil moisture data.'

        DO is = 1, nstyp1+1
          DO k = 1,nzsoil
            DO j = 1,nylg-1
              DO i = 1,nxlg-1
                var4d(i,j,k,is) = net4d(i,j,k,is)
              END DO
            END DO
          END DO
          CALL edgfill(var4d(:,:,:,is),1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nzsoil,1,nzsoil)
        END DO
      END IF
      CALL mpisplit4d(var4d,nxin,nyin,nzsoil,nstyp1+1,qsoil_in)
    END IF

    IF (wcanpin == 1) THEN
      IF (myproc == 0) THEN
        CALL netread3d(sd_id,0,0,"WETCANP",nxlg-1,nylg-1,nstyp1+1,net3d)
        WRITE (6, '(1x,a)') 'Read in the canopy water amount data.'

        DO is = 1, nstyp1+1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3d(i,j,is) = net3d(i,j,is)
            END DO
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nstyp1+1,1,nstyp1+1)
      END IF
      CALL mpisplit3d(var3d,nxin,nyin,nstyp1+1,wetcanp_in)
    END IF

    IF (snowdin == 1) THEN
      IF (myproc == 0) THEN
        CALL netread2d(sd_id,0,0,"SNOWDPTH",nxlg-1,nylg-1,net3d)
        WRITE (6, '(1x,a)') 'Read in the snow depth data.'
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = net3d(i,j,1)
          END DO
        END DO
        CALL edgfill(var3d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var3d,nxin,nyin,snowdpth)
    END IF

    IF (stypin == 1) THEN
      IF (myproc == 0) THEN
        CALL netread3di(sd_id,0,0,"SOILTYP",nxlg-1,nylg-1,nstyp1,net3di)
        WRITE (6, '(1x,a)') 'Read in soil type of soil data.'

        DO is = 1,nstyp1
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3di(i,j,is) = net3di(i,j,is)
            END DO
          END DO
        END DO
        CALL iedgfill(var3di,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nstyp1,1,nstyp1)
      END IF
      CALL mpisplit3di(var3di,nxin,nyin,nstyp1,soiltyp_in)
    END IF

    IF(myproc == 0) CALL netclose(sd_id)

    DEALLOCATE(net3d, net3di)
    DEALLOCATE(net4d)

  END IF

  IF (stypin == 0) THEN

    WRITE (6,'(1x,2a)') "READSPLITSOIL: WARNING, no check made for ",   &
         "consistency between soil types in surface and soil data sets."

    nstypin = MIN(nstyp1, nstyps)

    IF (zpsoilin == 1) zpsoil (:,:,:) = zpsoil_in(:,:,:)
    IF (tsoilin  == 1) tsoil  (:,:,:,0:nstypin) = tsoil_in(:,:,:,0:nstypin)
    IF (qsoilin  == 1) qsoil  (:,:,:,0:nstypin) = qsoil_in(:,:,:,0:nstypin)
    IF (wcanpin  == 1) wetcanp(:,:,0:nstypin)   = wetcanp_in(:,:,0:nstypin)

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstyp1,nstyp,tsoil,qsoil,wetcanp)

  ELSE

    IF (nstyp1 == 1) THEN
      tsoil_in  (:,:,:,1) = tsoil_in  (:,:,:,0)
      qsoil_in  (:,:,:,1) = qsoil_in  (:,:,:,0)
      wetcanp_in(:,:,1)   = wetcanp_in(:,:,0)
    ENDIF

    IF (soilmodel_forced == 0) THEN

      CALL remap_soil_vars(nx,ny,nzsoil,nstyp1,nstyp,                   &
                         tsoil_in,qsoil_in,wetcanp_in,soiltyp_in,       &
                         tsfcin,tsoilin,wsfcin,wdpin,qsoilin,wcanpin,   &
                         intver,                                        &
                         tsoil,qsoil,wetcanp,soiltyp)

    ENDIF

  ENDIF

  DEALLOCATE (tsoil_in,stat=istat)
  DEALLOCATE (qsoil_in,stat=istat)
  DEALLOCATE (wetcanp_in,stat=istat)
  DEALLOCATE (soiltyp_in,stat=istat)

  DEALLOCATE (itmp,stat=istat)
  DEALLOCATE (var3d,stat=istat)
  DEALLOCATE (var3di,stat=istat)
  DEALLOCATE (var4d,stat=istat)

  IF (intver == intver410) THEN

    DEALLOCATE (tem1)

    IF (tsfcin /= tsoilin .OR. wsfcin /= wdpin) THEN

      WRITE (6,'(1x,a,a,a/,a/,a/)')  'READSPLITSOIL: WARNING: ',        &
                   'The soilvar data is of version ', fmtver410,        &
                '. The inconsistency flag between tsfcin and tsoilin, ',&
                ' or between wsfin and wdpin, may cause some problems. '
    END IF

    tsoilin = max(tsfcin,tsoilin)
    qsoilin = max(wsfcin,wdpin)

  END IF

  GO TO 999

  998   WRITE (6,'(/a,i2/a)')                                           &
         '     Read error in surface data file '                        &
         //soilfile//' with the I/O unit ',flunit,                      &
         'The model will STOP in subroutine READSPLITSOIL.'

  CALL arpsstop("arpsstop called from READSPLITSOIL reading surface data",1)

  999 CONTINUE

  RETURN
END SUBROUTINE readsplitsoil
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTSFCDT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE wrtsfcdt( nx,ny,nstyps, sfcoutfl, dx,dy,                     &
                  mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,  &
                  stypout,vtypout,laiout,rfnsout,vegout,ndviout,        &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write a surface property data .
!
!  NOTE:
!    Changes made here should also be made in subroutine wrtjoinsfcdt
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  2/20/94
!
!  MODIFICATION HISTORY:
!
!  2000/03/29 (Gene Bassett)
!  Removed the globcst.inc include.
!
!  2000/03/29 (Gene Bassett)
!  Added HDF4 format.
!
!  2004/08/15 (Yunheng Wang)
!  Added NetCDF 3.0 format.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of model grid points in the x-dir. (east/west)
!  ny       Number of model grid points in the y-dir. (north/south)
!  sfcoutfl Name of the output surface property data file
!
!  soiltyp  Soil type in model domain
!  vegtyp   Vegetation type in model domain
!  lai      Leaf Area Index in model domain
!  roufns   Surface roughness
!  veg      Vegetation fraction
!  ndvi     NDVI
!
!  OUTPUT:
!
!  Output written to the surface data file, sfcdtfl:
!
!  nx       Number of model grid points in the x-direction
!  ny       Number of model grid points in the y-direction
!
!  mapproj  Type of map projection used to setup the analysis grid.
!  trulat1  1st real true latitude of map projection.
!  trulat2  2nd real true latitude of map projection.
!  trulon   Real true longitude of map projection.
!  sclfct   Map scale factor. At latitude = trulat1 and trulat2
!
!  dx       Model grid spacing in the x-direction east-west (meters)
!  dy       Model grid spacing in the y-direction east-west (meters)
!  ctrlat    Lat. at the origin of the model grid (deg. N)
!  ctrlon    Lon. at the origin of the model grid (deg. E)
!
!  stypout  Flag for output of soil type
!  vtypout  Flag for output of vegetation type
!  laiout   Flag for output of Leaf Area Index
!  rfnsout  Flag for output of surface roughness
!  vegout   Flag for output of vegetation fraction
!  ndviout  Flag for output of NDVI
!
!  soiltyp  Soil type in model domain
!  stypfrct Fraction of each soil type in each grid box
!  vegtyp   Vegetation type in model domain
!  lai      Leaf Area Index in model domain
!  roufns   Surface roughness
!  veg      Vegetation fraction
!  ndvi     NDVI
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx              ! Number of grid points in the x-direction
  INTEGER :: ny              ! Number of grid points in the y-direction
  INTEGER :: nstyps          ! Max number of soil types in a grid box
  CHARACTER (LEN=* ) :: sfcoutfl ! Surface property data file name

  REAL :: dx
  REAL :: dy
  INTEGER :: mapproj       ! Map projection
  REAL :: trulat1          ! 1st real true latitude of map projection
  REAL :: trulat2          ! 2nd real true latitude of map projection
  REAL :: trulon           ! Real true longitude of map projection
  REAL :: sclfct           ! Map scale factor
  REAL :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon           ! Center longitude of the model domain (deg. E)

  INTEGER :: stypout         ! Flag for output of soiltyp
  INTEGER :: vtypout         ! Flag for output of vegtyp
  INTEGER :: laiout          ! Flag for output of lai
  INTEGER :: rfnsout         ! Flag for output of roufns
  INTEGER :: vegout          ! Flag for output of veg
  INTEGER :: ndviout         ! Flag for output of ndvi

  INTEGER :: soiltyp (nx,ny,nstyps)  ! Soil type in model domain
  REAL    :: stypfrct(nx,ny,nstyps)  ! Fraction of soil types
  INTEGER :: vegtyp  (nx,ny)         ! Vegetation type in model domain

  REAL :: lai    (nx,ny)     ! Leaf Area Index in model domain
  REAL :: roufns (nx,ny)     ! Surface roughness
  REAL :: veg    (nx,ny)     ! Vegetation fraction
  REAL :: ndvi   (nx,ny)     ! NDVI
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: flunit
  INTEGER :: idummy
  REAL    :: rdummy
  INTEGER :: ierr

  INTEGER :: i,j,is

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
  REAL :: atmp1(1),atmp2(1)
         ! unused arrays in hdf routines since NO COMPRESSION

  INTEGER :: stat, sd_id

  INTEGER, ALLOCATABLE :: var3di(:,:,:)
  REAL,    ALLOCATABLE :: var3d (:,:,:)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  idummy = 0
  rdummy = 0.0

  IF (sfcdmp == 0) RETURN

  WRITE (6,'(1x,/,2a)') 'WRTSFCDT: Opening file ',trim(sfcoutfl)

!-----------------------------------------------------------------------
!
!  Write out in Fortran unformatted.
!
!-----------------------------------------------------------------------

  IF (sfcdmp == 1) THEN

    CALL getunit( flunit )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(sfcoutfl, '-F f77 -N ieee', ierr)

    OPEN (UNIT = flunit, FILE = trim(sfcoutfl),STATUS = 'unknown',    &
          FORM = 'unformatted', ACCESS = 'sequential')

    WRITE (flunit) nx, ny

    WRITE (flunit) mapproj, stypout, vtypout, laiout, rfnsout,        &
                   vegout,  ndviout, nstyp,   idummy, idummy,         &
                   idummy,  idummy,  idummy,  idummy, idummy,         &
                   idummy,  idummy,  idummy,  idummy, idummy

    WRITE (flunit) dx,      dy,     ctrlat, ctrlon, trulat1,          &
                   trulat2, trulon, sclfct, rdummy, rdummy,           &
                   rdummy,  rdummy, rdummy, rdummy, rdummy,           &
                   rdummy,  rdummy, rdummy, rdummy, rdummy

    IF ( stypout /= 0 ) THEN

      IF ( nstyp <= 1 ) THEN
        WRITE (flunit) ((soiltyp(i,j,1),i=1,nx),j=1,ny)
      ELSE
        DO is=1,nstyp
          WRITE (flunit) ((soiltyp (i,j,is),i=1,nx),j=1,ny)
          WRITE (flunit) ((stypfrct(i,j,is),i=1,nx),j=1,ny)
        END DO
      END IF

    END IF

    IF ( vtypout /= 0 ) WRITE (flunit) vegtyp
    IF ( laiout  /= 0 ) WRITE (flunit) lai
    IF ( rfnsout /= 0 ) WRITE (flunit) roufns
    IF ( vegout  /= 0 ) WRITE (flunit) veg
    IF ( ndviout /= 0 ) WRITE (flunit) ndvi

    CLOSE ( flunit )
    CALL retunit ( flunit )

  ELSE IF (sfcdmp == 3) THEN

!-----------------------------------------------------------------------
!
!  Write out in HDF4.
!
!-----------------------------------------------------------------------

    CALL hdfopen(trim(sfcoutfl), 2, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "WRTSFCDT: ERROR creating HDF4 file: ",             &
                  trim(sfcoutfl)
      CALL arpsstop("ERROR creating surface data file.",1)
    END IF

    ! Make sure nstyp & stypfrct are valid for one soil type
    IF (nstyp == 0) nstyp = 1
    IF (nstyp == 1) stypfrct = 1

    CALL hdfwrti(sd_id, 'nstyp', nstyp, stat)

    CALL hdfwrti(sd_id, 'nx', nx, stat)
    CALL hdfwrti(sd_id, 'ny', ny, stat)
    CALL hdfwrtr(sd_id, 'dx', dx, stat)
    CALL hdfwrtr(sd_id, 'dy', dy, stat)
    CALL hdfwrti(sd_id, 'mapproj', mapproj, stat)
    CALL hdfwrtr(sd_id, 'trulat1', trulat1, stat)
    CALL hdfwrtr(sd_id, 'trulat2', trulat2, stat)
    CALL hdfwrtr(sd_id, 'trulon', trulon, stat)
    CALL hdfwrtr(sd_id, 'sclfct', sclfct, stat)
    CALL hdfwrtr(sd_id, 'ctrlat', ctrlat, stat)
    CALL hdfwrtr(sd_id, 'ctrlon', ctrlon, stat)

    IF ( stypout /= 0 ) THEN

      CALL hdfwrt3di(soiltyp,nx,ny,nstyp,sd_id,0,0,                &
                     'soiltyp','Soil type','index')
      CALL hdfwrt3d(stypfrct,nx,ny,nstyp,sd_id,0,0,                &
                    'stypfrct','Soil type fractional coverage',    &
                    'fraction',itmp,atmp1,atmp2)
    END IF

    IF ( vtypout /= 0 ) THEN
      CALL hdfwrt2di(vegtyp,nx,ny,sd_id,0,0,                          &
                     'vegtyp','Vegetation type','index')
    END IF

    IF ( laiout  /= 0 ) THEN
      CALL hdfwrt2d(lai,nx,ny,sd_id,0,0,                              &
                    'lai','Leaf Area Index','index',itmp)
    END IF

    IF ( rfnsout /= 0 ) THEN
      CALL hdfwrt2d(roufns,nx,ny,sd_id,0,0,                           &
                    'roufns','Surface roughness','0-1',itmp)
    END IF

    IF ( vegout  /= 0 ) THEN
      CALL hdfwrt2d(veg,nx,ny,sd_id,0,0,                              &
                    'veg','Vegetation fraction','fraction',itmp)
    END IF

    IF ( ndviout /= 0 ) THEN
      CALL hdfwrt2d(ndvi,nx,ny,sd_id,0,0, 'ndvi',                     &
               'Normalized differential vegetation index','index',itmp)
    END IF

    CALL hdfclose(sd_id,stat)
    IF (stat /= 0) THEN
      WRITE (6,*) "WRTSFCDT: ERROR on closing file ",trim(sfcoutfl),    &
                  " (status",stat,")"
    END IF

  ELSE IF (sfcdmp == 7) THEN
!-----------------------------------------------------------------------
!
!  Write out in NetCDF format
!
!-----------------------------------------------------------------------

    ALLOCATE(var3di(nx-1,ny-1,nstyp), STAT = stat)
    ALLOCATE(var3d (nx-1,ny-1,nstyp), STAT = stat)

    ! Make sure nstyp & stypfrct are valid for one soil type
    IF (nstyp == 0) nstyp = 1
    IF (nstyp == 1) stypfrct = 1

!-----------------------------------------------------------------------
!
!  Define surface file dimension and variables
!
!-----------------------------------------------------------------------

    CALL netopen(TRIM(sfcoutfl), 'C' , sd_id)

    CALL net_define_sfc(sd_id, nx,ny,nstyp, dx,dy,mapproj,sclfct,       &
          trulat1,trulat2,trulon,ctrlat,ctrlon,                         &
          stypout,vtypout,laiout,rfnsout,vegout,ndviout,stat)

    IF ( stypout /= 0 ) THEN
      DO is = 1,nstyp
        DO j = 1,ny-1
          DO i = 1,nx-1
            var3di(i,j,is) = soiltyp(i,j,is)
            var3d (i,j,is) = stypfrct(i,j,is)
          END DO
        END DO
      END DO
      CALL netwrt3di(sd_id,0,0,'SOILTYP',var3di,nx-1,ny-1,nstyp)
      CALL netwrt3d (sd_id,0,0,'STYPFRCT',var3d,nx-1,ny-1,nstyp)
    END IF

    IF ( vtypout /= 0 ) THEN
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3di(i,j,1) = vegtyp(i,j)
        END DO
      END DO
      CALL netwrt2di(sd_id,0,0,'VEGTYP',var3di,nx-1,ny-1)
    END IF

    IF ( laiout  /= 0 ) THEN
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3d(i,j,1) = lai(i,j)
        END DO
      END DO
      CALL netwrt2d(sd_id,0,0,'LAI',var3d,nx-1,ny-1)
    END IF

    IF ( rfnsout /= 0 ) THEN
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3d(i,j,1) = roufns(i,j)
        END DO
      END DO
      CALL netwrt2d(sd_id,0,0,'ROUFNS',var3d,nx-1,ny-1)
    END IF

    IF ( vegout  /= 0 ) THEN
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3d(i,j,1) = veg(i,j)
        END DO
      END DO
      CALL netwrt2d(sd_id,0,0,'VEG',var3d,nx-1,ny-1)
    END IF

    IF ( ndviout /= 0 ) THEN
      DO j = 1,ny-1
        DO i = 1,nx-1
          var3d(i,j,1) = ndvi(i,j)
        END DO
      END DO
      CALL netwrt2d(sd_id,0,0,'NDVI',var3d,nx-1,ny-1)
    END IF

    CALL netclose(sd_id)

    DEALLOCATE(var3d,var3di)

  ELSE

    ! alternate dump format ...
    WRITE(6,*) 'The supported surface data dump format are ',           &
               'binary (sfcdmp=1) and HDF4 no compressed (sfcdmp = 3).'
    CALL arpsstop('Surface data dump format is not supported.',1)

  END IF

  RETURN
END SUBROUTINE wrtsfcdt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTJOINSFCDT               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtjoinsfcdt( nx,ny,nstyps, sfcoutfl, dx,dy,                 &
           mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,         &
           stypout,vtypout,laiout,rfnsout,vegout,ndviout,               &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Join and Write a surface property data .
!
!  Note: Changes made here must also be made in subroutine wrtsfcdt.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  8/20/03
!  Modified from subroutine wrtsfcdt.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of model grid points in the x-dir. (east/west)
!  ny       Number of model grid points in the y-dir. (north/south)
!  sfcoutfl Name of the output surface property data file
!
!  soiltyp  Soil type in model domain
!  vegtyp   Vegetation type in model domain
!  lai      Leaf Area Index in model domain
!  roufns   Surface roughness
!  veg      Vegetation fraction
!  ndvi     NDVI
!
!  OUTPUT:
!
!  Output written to the surface data file, sfcdtfl:
!
!  nx       Number of model grid points in the x-direction
!  ny       Number of model grid points in the y-direction
!
!  mapproj  Type of map projection used to setup the analysis grid.
!  trulat1  1st real true latitude of map projection.
!  trulat2  2nd real true latitude of map projection.
!  trulon   Real true longitude of map projection.
!  sclfct   Map scale factor. At latitude = trulat1 and trulat2
!
!  dx       Model grid spacing in the x-direction east-west (meters)
!  dy       Model grid spacing in the y-direction east-west (meters)
!  ctrlat   Lat. at the origin of the model grid (deg. N)
!  ctrlon   Lon. at the origin of the model grid (deg. E)
!
!  stypout  Flag for output of soil type
!  vtypout  Flag for output of vegetation type
!  laiout   Flag for output of Leaf Area Index
!  rfnsout  Flag for output of surface roughness
!  vegout   Flag for output of vegetation fraction
!  ndviout  Flag for output of NDVI
!
!  soiltyp  Soil type in model domain
!  stypfrct Fraction of each soil type in each grid box
!  vegtyp   Vegetation type in model domain
!  lai      Leaf Area Index in model domain
!  roufns   Surface roughness
!  veg      Vegetation fraction
!  ndvi     NDVI
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx              ! Number of grid points in the x-direction
  INTEGER :: ny              ! Number of grid points in the y-direction
  INTEGER :: nstyps          ! Max number of soil types in a grid box
  CHARACTER (LEN=* ) :: sfcoutfl ! Surface property data file name

  REAL :: dx
  REAL :: dy
  INTEGER :: mapproj       ! Map projection
  REAL :: trulat1          ! 1st real true latitude of map projection
  REAL :: trulat2          ! 2nd real true latitude of map projection
  REAL :: trulon           ! Real true longitude of map projection
  REAL :: sclfct           ! Map scale factor
  REAL :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon           ! Center longitude of the model domain (deg. E)

  INTEGER :: stypout         ! Flag for output of soiltyp
  INTEGER :: vtypout         ! Flag for output of vegtyp
  INTEGER :: laiout          ! Flag for output of lai
  INTEGER :: rfnsout         ! Flag for output of roufns
  INTEGER :: vegout          ! Flag for output of veg
  INTEGER :: ndviout         ! Flag for output of ndvi

  INTEGER :: soiltyp (nx,ny,nstyps)  ! Soil type in model domain
  REAL    :: stypfrct(nx,ny,nstyps)  ! Fraction of soil types
  INTEGER :: vegtyp  (nx,ny)         ! Vegetation type in model domain

  REAL :: lai    (nx,ny)     ! Leaf Area Index in model domain
  REAL :: roufns (nx,ny)     ! Surface roughness
  REAL :: veg    (nx,ny)     ! Vegetation fraction
  REAL :: ndvi   (nx,ny)     ! NDVI
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: flunit
  INTEGER :: idummy
  REAL :: rdummy
  INTEGER :: ierr

  INTEGER :: i,j,is

  INTEGER :: nxlg, nylg
  INTEGER :: istatus

  REAL,    ALLOCATABLE :: out2d(:,:), stypfrctlg(:,:,:)
  INTEGER, ALLOCATABLE :: vegtyplg(:,:), soiltyplg(:,:,:)

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
  REAL :: atmp1(1),atmp2(1)
         ! unused arrays in hdf routines since NO COMPRESSION

  INTEGER :: stat, sd_id

  INTEGER, ALLOCATABLE :: var3di(:,:,:)
  REAL,    ALLOCATABLE :: var3d(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

  idummy = 0
  rdummy = 0.0

  IF (sfcdmp == 0) RETURN

  ALLOCATE(soiltyplg (nxlg,nylg,nstyps), STAT = istatus)
  ALLOCATE(stypfrctlg(nxlg,nylg,nstyps), STAT = istatus)
  ALLOCATE(vegtyplg  (nxlg,nylg), STAT = istatus)
  ALLOCATE(out2d     (nxlg,nylg), STAT = istatus)

  IF( myproc == 0 ) THEN
    WRITE (6,'(a,a)') 'WRTJOINSFCDT: Opening file ',trim(sfcoutfl)
  END IF

!-----------------------------------------------------------------------
!
!  Write out in Fortran unformatted.
!
!-----------------------------------------------------------------------

  IF (sfcdmp == 1) THEN

    IF( myproc == 0 ) THEN
      CALL getunit( flunit )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(sfcoutfl, '-F f77 -N ieee', ierr)

      OPEN (UNIT = flunit, FILE = trim(sfcoutfl),STATUS = 'unknown',    &
            FORM = 'unformatted', ACCESS = 'sequential')

      WRITE (flunit) nxlg, nylg

      WRITE (flunit) mapproj, stypout, vtypout, laiout, rfnsout,        &
                     vegout,  ndviout, nstyp,   idummy, idummy,         &
                     idummy,  idummy,  idummy,  idummy, idummy,         &
                     idummy,  idummy,  idummy,  idummy, idummy

      WRITE (flunit) dx,      dy,     ctrlat, ctrlon, trulat1,          &
                     trulat2, trulon, sclfct, rdummy, rdummy,           &
                     rdummy,  rdummy, rdummy, rdummy, rdummy,           &
                     rdummy,  rdummy, rdummy, rdummy, rdummy
    END IF

    IF ( stypout /= 0 ) THEN
        CALL mpimerge3di(soiltyp, nx,ny,nstyps,soiltyplg)
        CALL mpimerge3d (stypfrct,nx,ny,nstyps,stypfrctlg)
        IF(myproc == 0 ) THEN       ! only root processor write
          IF ( nstyp <= 1 ) THEN
            WRITE (flunit) ((soiltyplg(i,j,1),i=1,nxlg),j=1,nylg)
          ELSE
            DO is=1,nstyp
              WRITE (flunit) ((soiltyplg (i,j,is),i=1,nxlg),j=1,nylg)
              WRITE (flunit) ((stypfrctlg(i,j,is),i=1,nxlg),j=1,nylg)
            END DO
          END IF
        END IF   ! myproc == 0

    END IF

    IF ( vtypout /= 0 ) THEN
      CALL mpimerge2di(vegtyp,nx,ny,vegtyplg)
      IF(myproc == 0) WRITE (flunit) vegtyplg
    END IF

    IF ( laiout  /= 0 ) THEN
      CALL mpimerge2d(lai,nx,ny,out2d)
      IF(myproc == 0) WRITE (flunit) out2d
    END IF

    IF ( rfnsout /= 0 ) THEN
      CALL mpimerge2d(roufns,nx,ny,out2d)
      IF(myproc == 0) WRITE (flunit) out2d
    END IF

    IF ( vegout  /= 0 ) THEN
      CALL mpimerge2d(veg,nx,ny,out2d)
      IF(myproc == 0) WRITE (flunit) out2d
    END IF

    IF ( ndviout /= 0 ) THEN
      CALL mpimerge2d(ndvi,nx,ny,out2d)
      IF(myproc == 0) WRITE (flunit) out2d
    END IF

    IF( myproc == 0 ) THEN
      CLOSE ( flunit )
      CALL retunit ( flunit )
    END IF

  ELSE IF (sfcdmp == 3) THEN

!-----------------------------------------------------------------------
!
!  Write out in HDF4.
!
!-----------------------------------------------------------------------

    IF(myproc == 0 ) THEN
      CALL hdfopen(trim(sfcoutfl), 2, sd_id)
      IF (sd_id < 0) THEN
        WRITE (6,*) "WRTJOINSFCDT: ERROR creating HDF4 file: ",             &
                    trim(sfcoutfl)
        CALL arpsstop("ERROR creating surface data file.",1)
      END IF
    END IF

    ! Make sure nstyp & stypfrct are valid for one soil type
    IF (nstyp == 0) nstyp = 1
    IF (nstyp == 1) stypfrct = 1

    IF(myproc == 0 ) THEN
      CALL hdfwrti(sd_id, 'nstyp', nstyp, stat)

      CALL hdfwrti(sd_id, 'nx', nxlg, stat)
      CALL hdfwrti(sd_id, 'ny', nylg, stat)
      CALL hdfwrtr(sd_id, 'dx', dx, stat)
      CALL hdfwrtr(sd_id, 'dy', dy, stat)
      CALL hdfwrti(sd_id, 'mapproj', mapproj, stat)
      CALL hdfwrtr(sd_id, 'trulat1', trulat1, stat)
      CALL hdfwrtr(sd_id, 'trulat2', trulat2, stat)
      CALL hdfwrtr(sd_id, 'trulon', trulon, stat)
      CALL hdfwrtr(sd_id, 'sclfct', sclfct, stat)
      CALL hdfwrtr(sd_id, 'ctrlat', ctrlat, stat)
      CALL hdfwrtr(sd_id, 'ctrlon', ctrlon, stat)
    END IF

    IF ( stypout /= 0 ) THEN
      CALL mpimerge3di(soiltyp,nx,ny,nstyp,soiltyplg)
      CALL mpimerge3d (stypfrct,nx,ny,nstyp,stypfrctlg)
      IF(myproc == 0) THEN
        CALL hdfwrt3di(soiltyplg,nxlg,nylg,nstyp,sd_id,0,0,        &
                       'soiltyp','Soil type','index')
        CALL hdfwrt3d(stypfrctlg,nxlg,nylg,nstyp,sd_id,0,0,        &
                      'stypfrct','Soil type fractional coverage',  &
                      'fraction',itmp,atmp1,atmp2)
      END IF

    END IF

    IF ( vtypout /= 0 ) THEN
      CALL mpimerge2di(vegtyp,nx,ny,vegtyplg)
      IF(myproc == 0) CALL hdfwrt2di(vegtyplg,nxlg,nylg,sd_id,0,0,    &
                               'vegtyp','Vegetation type','index')
    END IF

    IF ( laiout  /= 0 ) THEN
      CALL mpimerge2d(lai,nx,ny,out2d)
      IF(myproc == 0) CALL hdfwrt2d(out2d,nxlg,nylg,sd_id,0,0,        &
                               'lai','Leaf Area Index','index',itmp)
    END IF

    IF ( rfnsout /= 0 ) THEN
        CALL mpimerge2d(roufns,nx,ny,out2d)
        IF(myproc == 0) CALL hdfwrt2d(out2d,nxlg,nylg,sd_id,0,0,        &
                               'roufns','Surface roughness','0-1',itmp)
    END IF

    IF ( vegout  /= 0 ) THEN
        CALL mpimerge2d(veg,nx,ny,out2d)
        IF(myproc == 0) CALL hdfwrt2d(out2d,nxlg,nylg,sd_id,0,0,        &
                            'veg','Vegetation fraction','fraction',itmp)
    END IF

    IF ( ndviout /= 0 ) THEN
        CALL mpimerge2d(ndvi,nx,ny,out2d)
        IF(myproc == 0) CALL hdfwrt2d(out2d,nxlg,nylg,sd_id,0,0,'ndvi',&
                 'Normalized differential vegetation index','index',itmp)
    END IF

    IF(myproc == 0 ) THEN
      CALL hdfclose(sd_id,stat)
      IF (stat /= 0) THEN
        WRITE (6,*) "WRTJOINSFCDT: ERROR on closing file ",trim(sfcoutfl),    &
                    " (status",stat,")"
      END IF
    END IF

  ELSE IF (sfcdmp == 7) THEN
!-----------------------------------------------------------------------
!
!  Write out in NetCDF format
!
!-----------------------------------------------------------------------

    ALLOCATE(var3di(nxlg-1,nylg-1,nstyp), STAT = stat)
    ALLOCATE(var3d (nxlg-1,nylg-1,nstyp), STAT = stat)

    ! Make sure nstyp & stypfrct are valid for one soil type
    IF (nstyp == 0) nstyp = 1
    IF (nstyp == 1) stypfrct = 1

    IF (myproc == 0) THEN

!-----------------------------------------------------------------------
!
!  Define surface file dimension and variables
!
!-----------------------------------------------------------------------

      CALL netopen(TRIM(sfcoutfl), 'C' , sd_id)
      CALL net_define_sfc(sd_id,nxlg,nylg,nstyp, dx,dy,mapproj,sclfct,  &
                   trulat1,trulat2,trulon,ctrlat,ctrlon,                &
                   stypout,vtypout,laiout,rfnsout,vegout,ndviout,stat)
    END IF

    IF ( stypout /= 0 ) THEN
      CALL mpimerge3di(soiltyp,nx,ny,nstyp,soiltyplg)
      CALL mpimerge3d (stypfrct,nx,ny,nstyp,stypfrctlg)
      IF (myproc == 0) THEN
        DO is = 1,nstyp
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3di(i,j,is) = soiltyplg(i,j,is)
              var3d(i,j,is)  = stypfrctlg(i,j,is)
            END DO
          END DO
        END DO
        CALL netwrt3di(sd_id,0,0,'SOILTYP',var3di,nxlg-1,nylg-1,nstyp)
        CALL netwrt3d (sd_id,0,0,'STYPFRCT',var3d,nxlg-1,nylg-1,nstyp)
      END IF
    END IF

    IF ( vtypout /= 0 ) THEN
      CALL mpimerge2di(vegtyp,nx,ny,vegtyplg)
      IF(myproc == 0) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3di(i,j,1) = vegtyplg(i,j)
          END DO
        END DO
        CALL netwrt2di(sd_id,0,0,'VEGTYP',var3di,nxlg-1,nylg-1)
      END IF
    END IF

    IF ( laiout  /= 0 ) THEN
      CALL mpimerge2d(lai,nx,ny,out2d)
      IF(myproc == 0) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = out2d(i,j)
          END DO
        END DO
        CALL netwrt2d(sd_id,0,0,'LAI',var3d,nxlg-1,nylg-1)
      END IF
    END IF

    IF ( rfnsout /= 0 ) THEN
      CALL mpimerge2d(roufns,nx,ny,out2d)
      IF(myproc == 0) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = out2d(i,j)
          END DO
        END DO
        CALL netwrt2d(sd_id,0,0,'ROUFNS',var3d,nxlg-1,nylg-1)
      END IF
    END IF

    IF ( vegout  /= 0 ) THEN
      CALL mpimerge2d(veg,nx,ny,out2d)
      IF(myproc == 0) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = out2d(i,j)
          END DO
        END DO
        CALL netwrt2d(sd_id,0,0,'VEG',var3d,nxlg-1,nylg-1)
      END IF
    END IF

    IF ( ndviout /= 0 ) THEN
      CALL mpimerge2d(ndvi,nx,ny,out2d)
      IF(myproc == 0) THEN
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var3d(i,j,1) = out2d(i,j)
          END DO
        END DO
        CALL netwrt2d(sd_id,0,0,'NDVI',var3d,nxlg-1,nylg-1)
      END IF
    END IF

    IF(myproc == 0) CALL netclose(sd_id)

    DEALLOCATE(var3d,var3di)

  ELSE

    ! alternate dump format ...
    WRITE(6,*) 'The supported surface data dump format are ',           &
               'binary (sfcdmp=1) and HDF4 no compressed (sfcdmp = 3).'
    CALL arpsstop('Surface data dump format is not supported.',1)

  END IF

  DEALLOCATE(soiltyplg,  STAT = istatus)
  DEALLOCATE(stypfrctlg, STAT = istatus)
  DEALLOCATE(vegtyplg,   STAT = istatus)
  DEALLOCATE(out2d,      STAT = istatus)

  RETURN
END SUBROUTINE wrtjoinsfcdt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE READSFCDT                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE readsfcdt( nx,ny,nstyps,sfcfile,dx,dy,                       &
           mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,         &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the surface data sets from file sfcfile.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/3/94
!
!  MODIFICATION HISTORY:
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D array, veg(nx,ny), to the soil and vegetation data
!  set file.
!
!  2000/03/29 (Gene Bassett)
!  Removed the globcst.inc include.
!
!  2000/03/29 (Gene Bassett)
!  Added HDF4 format.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of model grid points in the x-dir. (east/west)
!  ny       Number of model grid points in the y-dir. (north/south)
!
!  mapproj  Type of map projection used to setup the analysis grid.
!  trulat1  1st real true latitude of map projection.
!  trulat2  2nd real true latitude of map projection.
!  trulon   Real true longitude of map projection.
!  sclfct   Map scale factor. At latitude = trulat1 and trulat2
!
!  dx       Model grid spacing in the x-direction east-west (meters)
!  dy       Model grid spacing in the y-direction east-west (meters)
!  ctrlat    Lat. at the origin of the model grid (deg. N)
!  ctrlon    Lon. at the origin of the model grid (deg. E)
!
!  OUTPUT:
!
!  soiltyp  Soil type in model domain
!  vegtyp   Vegetation type in model domain
!  lai      Leaf Area Index in model domain
!  roufns   Surface roughness
!  veg      Vegetation fraction
!  ndvi     NDVI
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx              ! Number of grid points in the x-direction
  INTEGER :: ny              ! Number of grid points in the y-direction
  INTEGER :: nstyps          ! Max number of soil types in a grid box

  CHARACTER (LEN=*) :: sfcfile ! Name of the surface data file

  REAL :: dx
  REAL :: dy
  INTEGER :: mapproj       ! Map projection
  REAL :: trulat1          ! 1st real true latitude of map projection
  REAL :: trulat2          ! 2nd real true latitude of map projection
  REAL :: trulon           ! Real true longitude of map projection
  REAL :: sclfct           ! Map scale factor
  REAL :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon           ! Center longitude of the model domain (deg. E)

  INTEGER :: soiltyp(nx,ny,nstyps)  ! Soil type in model domain
  REAL    :: stypfrct(nx,ny,nstyps) ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)  ! Vegetation type in model domain

  REAL :: lai    (nx,ny)     ! Leaf Area Index in model domain
  REAL :: roufns (nx,ny)     ! NDVI in model domain
  REAL :: veg    (nx,ny)     ! Vegetation fraction
  REAL :: ndvi   (nx,ny)     ! NDVI
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxin,nyin
  REAL :: dxin,dyin
  INTEGER :: mprojin
  INTEGER :: nstyp1, nstypin
  REAL :: trlat1in, trlat2in, trlonin, sclfctin
  REAL :: ctrlonin, ctrlatin
  INTEGER :: stypin,vtypin,laiin,roufin,vegin,ndviin

  INTEGER :: idummy
  REAL :: rdummy

  INTEGER :: istat, ierr,i,j,is

  INTEGER :: ireturn

  CHARACTER (LEN=256) :: savename
  CHARACTER :: ayear*4, amm*2, aday*2

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
  REAL                               :: atmp1(1),atmp2(1)
         ! unused arrays in hdf routines since NO COMPRESSION

  INTEGER :: stat, sd_id

  INTEGER, ALLOCATABLE :: temi(:,:,:)
  REAL,    ALLOCATABLE :: temr(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Open the surface data file. Read the parameters first, check if
!  the data are consistant with the model. If everything is OK, then
!  read the surface data, soiltyp, vegtyp, lai, roufns, veg, and
!  ndvi.
!
!-----------------------------------------------------------------------
!

  IF (mp_opt > 0) THEN
    savename(1:256) = sfcfile(1:256)
    CALL gtsplitfn(savename,1,1,loc_x,loc_y,1,1,0,0,1,lvldbg,sfcfile,ireturn)
  END IF

  WRITE (6,'(1x,3a)') 'READSFCDT: reading in external supplied surface',&
                      'data from file - ',trim(sfcfile)

!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------

  IF (sfcfmt == 1) THEN

!-----------------------------------------------------------------------
!
!  Fortran unformatted dump.
!
!-----------------------------------------------------------------------

    CALL getunit( sfcunit )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(sfcfile, '-F f77 -N ieee', ierr)

    OPEN(UNIT=sfcunit,FILE=trim(sfcfile),FORM='unformatted',            &
         STATUS='old',IOSTAT=istat)

    IF( istat /= 0 ) THEN

      WRITE(6,'(/1x,a,i2,/1x,a/)')                                      &
          'Error occured when opening the surface data file '           &
          //sfcfile//' using FORTRAN unit ',sfcunit,                    &
          ' Program stopped in READSFCDT.'
      CALL arpsstop("arpsstop called from READSFCDT opening file",1)

    END IF

    WRITE(6,'(/1x,a,/1x,a,i2/)')                                        &
        'This run will start from an external supplied surface ',       &
        'data file '//sfcfile//' using FORTRAN unit ',sfcunit

    READ (sfcunit,ERR=998) nxin,nyin

    READ (sfcunit,ERR=998) mprojin,stypin,vtypin,laiin,roufin,          &
                         vegin,  ndviin,nstyp1,idummy,idummy,           &
                         idummy, idummy,idummy,idummy,idummy,           &
                         idummy, idummy,idummy,idummy,idummy

    READ (sfcunit,ERR=998) dxin,dyin, ctrlatin,ctrlonin,trlat1in,       &
                         trlat2in,trlonin,sclfctin,rdummy,rdummy,       &
                         rdummy,rdummy,rdummy,rdummy,rdummy,            &
                         rdummy,rdummy,rdummy,rdummy,rdummy

  ELSE IF (sfcfmt == 3) THEN

!-----------------------------------------------------------------------
!
!  HDF4 format.
!
!-----------------------------------------------------------------------

    CALL hdfopen(trim(sfcfile), 1, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "READSFCDT: ERROR opening ",                          &
                 trim(sfcfile)," for reading."
      GO TO 998
    END IF

    CALL hdfrdi(sd_id,"nstyp",nstyp1,istat)

    CALL hdfrdi(sd_id,"nx",nxin,istat)
    CALL hdfrdi(sd_id,"ny",nyin,istat)
    CALL hdfrdr(sd_id,"dx",dxin,istat)
    CALL hdfrdr(sd_id,"dy",dyin,istat)
    CALL hdfrdi(sd_id,"mapproj",mprojin,istat)
    CALL hdfrdr(sd_id,"trulat1",trlat1in,istat)
    CALL hdfrdr(sd_id,"trulat2",trlat2in,istat)
    CALL hdfrdr(sd_id,"trulon",trlonin,istat)
    CALL hdfrdr(sd_id,"sclfct",sclfctin,istat)
    CALL hdfrdr(sd_id,"ctrlat",ctrlatin,istat)
    CALL hdfrdr(sd_id,"ctrlon",ctrlonin,istat)

  ELSE IF (sfcfmt == 7) THEN

!-----------------------------------------------------------------------
!
!  NetCDF format.
!
!-----------------------------------------------------------------------

    CALL netopen(TRIM(sfcfile), 'R', sd_id)

    CALL net_get_sfc(sd_id,nxin,nyin,nstyp1,dxin,dyin,                  &
           mprojin,sclfctin,trlat1in,trlat2in,trlonin,ctrlatin,ctrlonin,&
           stypin,vtypin,laiin,roufin,vegin,ndviin,istat)

  ELSE IF (sfcfmt == 2) THEN    ! Direct output        (JAB)

    CALL initztime(ayear,amm,aday)

    CALL readjveg(nx,ny,nstyp,ayear,amm,aday,ztime,                     &
             soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi)

  ELSE

    ! alternate dump format ...

    WRITE(6,*) 'The supported surface data format are ',                &
               'binary (sfcfmt=1) and HDF4 no compressed (sfcfmt = 3).'
    CALL arpsstop('Surface data format is not supported.',1)

  END IF                     !sfcfmt loop

  nstyp1 = MAX( nstyp1, 1 )

!-----------------------------------------------------------------------
!
!  Check the data file for consistent grid parameters.
!
!-----------------------------------------------------------------------

  IF (sfcfmt /= 2) THEN

    CALL checkgrid2d(nx,ny,nxin,nyin,                                   &
                dx,dy,ctrlat,ctrlon,                                    &
                mapproj,trulat1,trulat2,trulon,sclfct,                  &
                dxin,dyin,ctrlatin,ctrlonin,                            &
                mprojin,trlat1in,trlat2in,trlonin,sclfctin,ireturn)

  IF (ireturn /= 0) THEN
    WRITE (6,*) "READSFCDT: ERROR, grid parameter mismatch"
    CALL arpsstop("arpsstop called from READSFCDT grid param. mismatch",1)
  END IF

  IF (myproc == 0) WRITE (6,'(a,a//a,I2/6(a,e15.8/))')                  &
      ' The map projection and griding information for the ',           &
      ' surface data: ',                                                &
      ' Projection:                 ', mprojin,                         &
      ' The 1st real true latitude: ', trlat1in,                        &
      ' The 2nd real true latitude: ', trlat2in,                        &
      ' The real true longitude:    ', trlonin,                         &
      ' Map scale factor:           ', sclfctin,                        &
      ' Latitude  at the origin:    ', ctrlatin,                        &
      ' Longitude at the origin:    ', ctrlonin


  ENDIF

!-----------------------------------------------------------------------
!
!  Read in the surface data from the surface data file.
!
!-----------------------------------------------------------------------

  IF (sfcfmt == 1) THEN    ! Fortran unformatted

    WRITE (6, '(a/a,i2/a,i2/a,i2/a,i2/a,i2/a,i2)')                        &
      ' Surface data flags for: ',                                      &
      '        soiltyp --', stypin,                                     &
      '         vegtyp --', vtypin,                                     &
      '            lai --', laiin,                                      &
      '         roufns --', roufin,                                     &
      '            veg --', vegin,                                      &
      '           ndvi --', ndviin

    WRITE (6, '(a/a,i2)')                                                 &
      ' Number of soil types in each grid box:',                        &
      '          nstyp --', nstyp

    IF(stypin == 1) THEN
      IF ( nstyp1 == 1 ) THEN
        WRITE (6, '(a)') 'Read in the soil type data'
        READ (sfcunit,ERR=998) ((soiltyp(i,j,1),i=1,nx),j=1,ny)
        DO j=1,ny
          DO i=1,nx
            stypfrct(i,j,1) = 1.0
          END DO
        END DO
      ELSE
        DO is=1,nstyp1
          IF (is <= nstyps) THEN
            WRITE (6, '(a)') 'Read in the soil type data'
            READ (sfcunit,ERR=998) ((soiltyp(i,j,is),i=1,nx),j=1,ny)
            WRITE (6, '(a)') 'Read in the fraction of soil types'
            READ (sfcunit,ERR=998) ((stypfrct(i,j,is),i=1,nx),j=1,ny)
          ELSE
            READ (sfcunit,ERR=998)
            READ (sfcunit,ERR=998)
          ENDIF
        END DO
      END IF
    END IF

    IF(vtypin == 1) THEN
      WRITE (6, '(a)') 'Read in the vegetation type data'
      READ (sfcunit,ERR=998) vegtyp
    END IF

    IF(laiin == 1) THEN
      WRITE (6, '(a)') 'Read in the Leaf Area Index data'
      READ (sfcunit,ERR=998) lai
    END IF

    IF(roufin == 1) THEN
      WRITE (6, '(a)') 'Read in the surface roughness data'
      READ (sfcunit,ERR=998) roufns
    END IF

    IF(vegin == 1) THEN
      WRITE (6, '(a)') 'Read in the vegetatin fraction data'
      READ (sfcunit,ERR=998) veg
    END IF

    IF (ndviin == 1) THEN
      WRITE (6, '(a)') 'Read in the NDVI data'
      READ (sfcunit,ERR=998) ndvi
    END IF

    CLOSE ( sfcunit )
    CALL retunit ( sfcunit )

  ELSE IF (sfcfmt == 3) THEN

    nstypin = MIN(nstyp1, nstyps)

    CALL hdfrd3di(sd_id,"soiltyp",nx,ny,nstypin,soiltyp,stat)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(1x,a)') 'Read in soiltyp'
      stypin = 1
    ELSE
      WRITE (6, '(1x,a)') 'Variable soiltyp is not in the data set.'
      stypin = 0
    END IF
    CALL hdfrd3d(sd_id,"stypfrct",nx,ny,nstypin,                    &
        stypfrct,stat,itmp,atmp1,atmp2)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(1x,a)') 'Read in stypfrct'
    ELSE
      WRITE (6, '(1x,a)') 'Variable stypfrct is not in the data set.'
      stypfrct(:,:,1) = 1.
      stypfrct(:,:,2:nstypin) = 0.
    END IF

    CALL hdfrd2di(sd_id,"vegtyp",nx,ny,vegtyp,stat)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(1x,a)') 'Read in vegtyp'
      vtypin = 1
    ELSE
      WRITE (6, '(1x,a)') 'Variable vegtyp is not in the data set.'
      vtypin = 0
    END IF

    CALL hdfrd2d(sd_id,"lai",nx,ny,lai,stat,itmp)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(1x,a)') 'Read in lai'
      laiin = 1
    ELSE
      WRITE (6, '(1x,a)') 'Variable lai is not in the data set.'
      laiin = 0
    END IF

    CALL hdfrd2d(sd_id,"roufns",nx,ny,roufns,stat,itmp)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(1x,a)') 'Read in roufns'
      roufin = 1
    ELSE
      WRITE (6, '(1x,a)') 'Variable roufns is not in the data set.'
      roufin = 0
    END IF

    CALL hdfrd2d(sd_id,"veg",nx,ny,veg,stat,itmp)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(1x,a)') 'Read in veg'
      vegin = 1
    ELSE
      WRITE (6, '(1x,a)') 'Variable veg is not in the data set.'
      vegin = 0
    END IF

    CALL hdfrd2d(sd_id,"ndvi",nx,ny,ndvi,stat,itmp)
    IF (stat > 1) GO TO 998
    IF (stat == 0) THEN
      WRITE (6, '(1x,a)') 'Read in ndvi'
      ndviin = 1
    ELSE
      WRITE (6, '(1x,a)') 'Variable ndvi is not in the data set.'
      ndviin = 0
    END IF

    CALL hdfclose(sd_id, stat)

  ELSE IF (sfcfmt == 7) THEN    ! NetCDF format

    WRITE (6, '(a/a,i2/a,i2/a,i2/a,i2/a,i2/a,i2)')                      &
                            ' Surface data flags for: ',                &
                            '        soiltyp --', stypin,               &
                            '         vegtyp --', vtypin,               &
                            '            lai --', laiin,                &
                            '         roufns --', roufin,               &
                            '            veg --', vegin,                &
                            '           ndvi --', ndviin

    WRITE (6, '(a/a,i2)') ' Number of soil types in each grid box:',    &
                          '          nstyp --', nstyp

    nstypin = MIN(nstyp1,nstyps)

    ALLOCATE(temr(nx-1,ny-1,nstyp1), STAT = stat)
    ALLOCATE(temi(nx-1,ny-1,nstyp1), STAT = stat)

    IF(stypin == 1) THEN
      CALL netread3di(sd_id,0,0,'SOILTYP',nx-1,ny-1,nstyp1,temi)
      CALL netread3d (sd_id,0,0,'STYPFRCT',nx-1,ny-1,nstyp1,temr)
      DO is = 1, nstypin
        DO j = 1, ny-1
          DO i = 1, nx-1
            soiltyp (i,j,is) = temi(i,j,is)
            stypfrct(i,j,is) = temr(i,j,is)
          END DO
        END DO
      END DO
      CALL iedgfill(soiltyp, 1,nx,1,nx-1,1,ny,1,ny-1,1,nstypin,1,nstypin)
      CALL  edgfill(stypfrct,1,nx,1,nx-1,1,ny,1,ny-1,1,nstypin,1,nstypin)
    END IF

    IF(vtypin == 1) THEN
      CALL netread2di(sd_id,0,0,'VEGTYP',nx-1,ny-1,temi)
      DO j = 1, ny-1
        DO i = 1, nx-1
          vegtyp (i,j) = temi(i,j,1)
        END DO
      END DO
      CALL iedgfill(vegtyp, 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    IF(laiin == 1) THEN
      CALL netread2d(sd_id,0,0,'LAI',nx-1,ny-1,temr)
      DO j = 1, ny-1
        DO i = 1, nx-1
          lai (i,j) = temr(i,j,1)
        END DO
      END DO
      CALL edgfill(lai,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    IF(roufin == 1) THEN
      CALL netread2d(sd_id,0,0,'ROUFNS',nx-1,ny-1,temr)
      DO j = 1, ny-1
        DO i = 1, nx-1
          roufns (i,j) = temr(i,j,1)
        END DO
      END DO
      CALL edgfill(roufns,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    IF(vegin == 1) THEN
      CALL netread2d(sd_id,0,0,'VEG',nx-1,ny-1,temr)
      DO j = 1, ny-1
        DO i = 1, nx-1
          veg (i,j) = temr(i,j,1)
        END DO
      END DO
      CALL edgfill(veg,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    IF (ndviin == 1) THEN
      CALL netread2d(sd_id,0,0,'NDVI',nx-1,ny-1,temr)
      DO j = 1, ny-1
        DO i = 1, nx-1
          ndvi (i,j) = temr(i,j,1)
        END DO
      END DO
      CALL edgfill(ndvi,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END IF

    CALL netclose(sd_id)

    DEALLOCATE(temi, temr)

    ! alternate dump format ...

  ELSE IF (sfcfmt == 2) THEN    !Direct output from OASIS data

    CALL initztime(ayear,amm,aday)

    CALL readjveg(nx,ny,nstyp,ayear,amm,aday,ztime,       &
             soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi)

  END IF

  CALL fix_stypfrct_nstyp(nx,ny,nstyp1,nstyp,stypfrct)

  IF(stypin == 0) THEN
    WRITE(6,'(/1x,a,a)')                                                &
        'Soil type data not present in the surface propery ',           &
        'data file, set it to styp.'
    DO i=1,nx
      DO j=1,ny
        soiltyp(i,j,1) = styp
      END DO
    END DO
  END IF

  IF(vtypin == 0) THEN
    WRITE(6,'(/1x,a,a)')                                                &
        'Vegetation type data not present in the surface propery ',     &
        'data file, set it to vtyp.'

    DO i=1,nx-1
      DO j=1,ny-1
        vegtyp(i,j) = vtyp
      END DO
    END DO
  END IF

  IF(laiin == 0) THEN
    WRITE(6,'(/1x,a,a)')                                                &
        'Leaf Area Index data not present in the surface propery ',     &
        'data file, set it to lai0.'

    DO i=1,nx-1
      DO j=1,ny-1
        lai(i,j) = lai0
      END DO
    END DO
  END IF

  IF(roufin == 0) THEN
    WRITE(6,'(/1x,a,a)')                                                &
        'Roughness data not present in the surface propery ',           &
        'data file, set it to roufns0.'

    DO i=1,nx-1
      DO j=1,ny-1
        roufns(i,j) = roufns0
      END DO
    END DO
  END IF

  IF(vegin == 0) THEN
    WRITE(6,'(/1x,a,a)')                                                &
        'Vegetation fraction not present in the surface propery ',      &
        'data file, set it to veg0.'

    DO i=1,nx-1
      DO j=1,ny-1
        veg(i,j) = veg0
      END DO
    END DO
  END IF

  IF (mp_opt > 0) sfcfile(1:256) = savename(1:256)

  RETURN

  998   WRITE (6,'(/a,i2/a)')                                           &
         'READSFCDT: Read error in surface data file '                  &
         //sfcfile//' with the I/O unit ',sfcunit,                      &
         'The model will STOP in subroutine READSFCDT.'

  CALL arpsstop("arpsstop called from READSFCDT reading sfc file",1)

END SUBROUTINE readsfcdt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE READSPLITSFCDT            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readsplitsfcdt( nx,ny,nstyps,sfcfile,dx,dy,                  &
           mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,         &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the surface data sets from file sfcfile. Split and scatter the
!  data to MP processors.
!
!  NOTE:
!
!  Changes in one of the subroutine readsfcdt or readsplitsfcdt should
!  also be maken in the other.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  8/30/2002
!  Based on subroutine readsfcdt.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of model grid points in the x-dir. (east/west)
!  ny       Number of model grid points in the y-dir. (north/south)
!
!  mapproj  Type of map projection used to setup the analysis grid.
!  trulat1  1st real true latitude of map projection.
!  trulat2  2nd real true latitude of map projection.
!  trulon   Real true longitude of map projection.
!  sclfct   Map scale factor. At latitude = trulat1 and trulat2
!
!  dx       Model grid spacing in the x-direction east-west (meters)
!  dy       Model grid spacing in the y-direction east-west (meters)
!  ctrlat    Lat. at the origin of the model grid (deg. N)
!  ctrlon    Lon. at the origin of the model grid (deg. E)
!
!  OUTPUT:
!
!  soiltyp  Soil type in model domain
!  vegtyp   Vegetation type in model domain
!  lai      Leaf Area Index in model domain
!  roufns   Surface roughness
!  veg      Vegetation fraction
!  ndvi     NDVI
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx              ! Number of grid points in the x-direction
  INTEGER :: ny              ! Number of grid points in the y-direction
  INTEGER :: nstyps          ! Max number of soil types in a grid box

  CHARACTER (LEN=*) :: sfcfile ! Name of the surface data file

  REAL :: dx
  REAL :: dy
  INTEGER :: mapproj       ! Map projection
  REAL :: trulat1          ! 1st real true latitude of map projection
  REAL :: trulat2          ! 2nd real true latitude of map projection
  REAL :: trulon           ! Real true longitude of map projection
  REAL :: sclfct           ! Map scale factor
  REAL :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon           ! Center longitude of the model domain (deg. E)

  INTEGER :: soiltyp(nx,ny,nstyps)  ! Soil type in model domain
  REAL :: stypfrct(nx,ny,nstyps)    ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)  ! Vegetation type in model domain

  REAL :: lai    (nx,ny)     ! Leaf Area Index in model domain
  REAL :: roufns (nx,ny)     ! NDVI in model domain
  REAL :: veg    (nx,ny)     ! Vegetation fraction
  REAL :: ndvi   (nx,ny)     ! NDVI
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxin,nyin
  REAL    :: dxin,dyin
  INTEGER :: mprojin
  INTEGER :: nstyp1, nstypin
  REAL    :: trlat1in, trlat2in, trlonin, sclfctin
  REAL    :: ctrlonin, ctrlatin
  INTEGER :: stypin,vtypin,laiin,roufin,vegin,ndviin

  REAL*4  :: dxin4,dyin4
  REAL*4  :: trlat1in4, trlat2in4, trlonin4, sclfctin4
  REAL*4  :: ctrlonin4, ctrlatin4
  REAL*4  :: rdummy4

  INTEGER :: idummy
  REAL    :: rdummy

  INTEGER :: istat, ierr,i,j,is

  INTEGER :: ireturn

  CHARACTER :: ayear*4, amm*2, aday*2

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
  REAL                               :: atmp1(1),atmp2(1)
         ! unused arrays in hdf routines since NO COMPRESSION

  INTEGER :: stat, sd_id

  INTEGER :: nxlg, nylg
  INTEGER, ALLOCATABLE :: var2di(:,:), var3di(:,:,:)
  REAL,    ALLOCATABLE :: var2d(:,:), var3d(:,:,:)

  REAL,    ALLOCATABLE :: temr(:,:,:)
  INTEGER, ALLOCATABLE :: temi(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3

!-----------------------------------------------------------------------
!
!  Open the surface data file. Read the parameters first, check if
!  the data are consistant with the model. If everything is OK, then
!  read the surface data, soiltyp, vegtyp, lai, roufns, veg, and
!  ndvi.
!
!-----------------------------------------------------------------------
!

  IF (myproc == 0) THEN

    WRITE (6,'(1x,3a)') 'READSPLITSFCDT: reading in external supplied surface ',  &
                'data from file - ',trim(sfcfile)

!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------

    IF (sfcfmt == 1) THEN

!-----------------------------------------------------------------------
!
!  Fortran unformatted dump.
!
!-----------------------------------------------------------------------

      CALL getunit( sfcunit )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(sfcfile, '-F f77 -N ieee', ierr)

      OPEN(UNIT=sfcunit,FILE=trim(sfcfile),FORM='unformatted',            &
           STATUS='old',IOSTAT=istat)

      IF( istat /= 0 ) THEN

        WRITE(6,'(/1x,a,i2,/1x,a/)')                                      &
            'Error occured when opening the surface data file '           &
            //sfcfile//' using FORTRAN unit ',sfcunit,                    &
            ' Program stopped in READSPLITSFCDT.'
        CALL arpsstop("arpsstop called from READSPLITSFCDT opening file",1)

      END IF

      WRITE(6,'(/1x,a,/1x,a,i2/)')                                        &
          'This run will start from an external supplied surface ',       &
          'data file '//sfcfile//' using FORTRAN unit ',sfcunit

      READ (sfcunit,ERR=998) nxin,nyin

      READ (sfcunit,ERR=998) mprojin, stypin, vtypin,  laiin, roufin,   &
                               vegin, ndviin, nstyp1, idummy, idummy,   &
                              idummy, idummy, idummy, idummy, idummy,   &
                              idummy, idummy, idummy, idummy, idummy

      READ (sfcunit,ERR=998) dxin4,   dyin4, ctrlatin4,ctrlonin4,trlat1in4,&
                         trlat2in4,trlonin4, sclfctin4,  rdummy4,  rdummy4,&
                           rdummy4, rdummy4,   rdummy4,  rdummy4,  rdummy4,&
                           rdummy4, rdummy4,   rdummy4,  rdummy4,  rdummy4

       dxin = dxin4
       dyin = dyin4
       trlat1in = trlat1in4
       trlat2in = trlat2in4
       trlonin  = trlonin4
       sclfctin = sclfctin4
       ctrlonin = ctrlonin4
       ctrlatin = ctrlatin4

    ELSE IF (sfcfmt == 3) THEN  !  HDF4 format.

      CALL hdfopen(trim(sfcfile), 1, sd_id)
      IF (sd_id < 0) THEN
        WRITE (6,*) "READSPLITSFCDT: ERROR opening ",                &
                   trim(sfcfile)," for reading."
        GO TO 998
      END IF

      CALL hdfrdi(sd_id,"nstyp",nstyp1,istat)

      CALL hdfrdi(sd_id,"nx",nxin,istat)
      CALL hdfrdi(sd_id,"ny",nyin,istat)
      CALL hdfrdr(sd_id,"dx",dxin,istat)
      CALL hdfrdr(sd_id,"dy",dyin,istat)
      CALL hdfrdi(sd_id,"mapproj",mprojin,istat)
      CALL hdfrdr(sd_id,"trulat1",trlat1in,istat)
      CALL hdfrdr(sd_id,"trulat2",trlat2in,istat)
      CALL hdfrdr(sd_id,"trulon",  trlonin,istat)
      CALL hdfrdr(sd_id,"sclfct", sclfctin,istat)
      CALL hdfrdr(sd_id,"ctrlat", ctrlatin,istat)
      CALL hdfrdr(sd_id,"ctrlon", ctrlonin,istat)

    ELSE IF (sfcfmt == 7) THEN

!-----------------------------------------------------------------------
!
!  NetCDF format.
!
!-----------------------------------------------------------------------

      CALL netopen(TRIM(sfcfile), 'R', sd_id)

      CALL net_get_sfc(sd_id,nxin,nyin,nstyp1,dxin,dyin,                &
          mprojin,sclfctin,trlat1in,trlat2in,trlonin,ctrlatin,ctrlonin, &
          stypin,vtypin,laiin,roufin,vegin,ndviin,istat)

    ELSE                      ! alternate dump format ...

      WRITE(6,*) 'READSPLITSFCDT only support binary and HDF',      &
                 'format right now. PROGRAM stopping!!!!'
      CALL arpsstop('Surface data format not supported',1)

    END IF                    !sfcfmt loop

    nxin = (nxin-3)/nproc_x+3
    nyin = (nyin-3)/nproc_y+3
    nstyp1 = MAX( nstyp1, 1 )

  END IF          ! myproc == 0

  IF (sfcfmt == 1 .OR. sfcfmt == 7) THEN
    CALL mpupdatei(stypin,1)
    CALL mpupdatei(vtypin,1)
    CALL mpupdatei(laiin,1)
    CALL mpupdatei(roufin,1)
    CALL mpupdatei(vegin,1)
    CALL mpupdatei(ndviin,1)
  END IF
  CALL mpupdatei(nstyp1,1)
  CALL mpupdatei(nxin,1)
  CALL mpupdatei(nyin,1)
  CALL mpupdater(dxin,1)
  CALL mpupdater(dyin,1)
  CALL mpupdatei(mprojin,1)
  CALL mpupdater(trlat1in,1)
  CALL mpupdater(trlat2in,1)
  CALL mpupdater(trlonin,1)
  CALL mpupdater(sclfctin,1)
  CALL mpupdater(ctrlatin,1)
  CALL mpupdater(ctrlonin,1)


!-----------------------------------------------------------------------
!
!  Check the data file for consistent grid parameters.
!
!-----------------------------------------------------------------------

  CALL checkgrid2d(nx,ny,nxin,nyin,                                     &
        dx,dy,ctrlat,ctrlon,mapproj,trulat1,trulat2,trulon,sclfct,      &
        dxin,dyin,ctrlatin,ctrlonin,mprojin,trlat1in,trlat2in,trlonin,sclfctin, &
        ireturn)

  IF (ireturn /= 0) THEN
    WRITE (6,*) "READSPLITSFCDT: ERROR, grid parameter mismatch"
    CALL arpsstop("arpsstop called from READSPLITSFCDT grid param. mismatch",1)
  END IF

  IF (myproc == 0) &

    WRITE (6,'(a,a//a,i2/6(a,e15.8/))')                                   &
        ' The map projection and griding information for the ',           &
        ' surface data: ',                                                &
        ' Projection:                 ', mprojin,                         &
        ' The 1st real true latitude: ', trlat1in,                        &
        ' The 2nd real true latitude: ', trlat2in,                        &
        ' The real true longitude:    ', trlonin,                         &
        ' Map scale factor:           ', sclfctin,                        &
        ' Latitude  at the origin:    ', ctrlatin,                        &
        ' Longitude at the origin:    ', ctrlonin


!-----------------------------------------------------------------------
!
!  Read in the surface data from the surface data file.
!
!-----------------------------------------------------------------------

  IF (sfcfmt == 1) THEN    ! Fortran unformatted

    ALLOCATE(var2di(nxlg,nylg))
    ALLOCATE(var2d(nxlg,nylg))
    ALLOCATE(var3di(nxlg,nylg,nstyps))
    ALLOCATE(var3d(nxlg,nylg,nstyps))

    IF (myproc == 0)  THEN

      WRITE (6, '(a/a,i2/a,i2/a,i2/a,i2/a,i2/a,i2)')                      &
        ' Surface data flags for: ',                                      &
        '        soiltyp --', stypin,                                     &
        '         vegtyp --', vtypin,                                     &
        '            lai --', laiin,                                      &
        '         roufns --', roufin,                                     &
        '            veg --', vegin,                                      &
        '           ndvi --', ndviin

      WRITE (6, '(a/a,i2)')                                               &
        ' Number of soil types in each grid box:',                        &
        '          nstyp --', nstyp

    END IF

    IF(stypin == 1) THEN
      IF ( nstyp1 == 1 ) THEN
        IF (myproc == 0) THEN
          WRITE (6, '(a)') 'Read in the soil type data'
          READ (sfcunit,ERR=998) ((var2di(i,j),i=1,nxlg),j=1,nylg)
        END IF
        CALL mpisplit3di(var2di,nx,ny,1,soiltyp(1,1,1))
        DO j=1,ny
          DO i=1,nx
            stypfrct(i,j,1) = 1.0
          END DO
        END DO
      ELSE
        DO is=1,nstyp1
          IF (is <= nstyps) THEN
            IF (myproc == 0) THEN
              WRITE (6, '(a)') 'Read in the soil type data'
              READ (sfcunit,ERR=998) ((var2di(i,j),i=1,nxlg),j=1,nylg)
              WRITE (6, '(a)') 'Read in the fraction of soil types'
              READ (sfcunit,ERR=998) ((var2d(i,j),i=1,nxlg),j=1,nylg)
            END IF
            CALL mpisplit3di(var2di,nx,ny,1,soiltyp(1,1,is))
            CALL mpisplit3d(var2d,nx,ny,1,stypfrct(1,1,is))
          ELSE
            IF (myproc == 0) THEN
              READ (sfcunit,ERR=998)
              READ (sfcunit,ERR=998)
            END IF  ! myproc == 0
          ENDIF
        END DO
      END IF
    END IF

    IF(vtypin == 1) THEN
      IF (myproc == 0) THEN
        WRITE (6, '(a)') 'Read in the vegetation type data'
        READ (sfcunit,ERR=998) var2di
      END IF
      CALL mpisplit3di(var2di,nx,ny,1,vegtyp)
    END IF

    IF(laiin == 1) THEN
      IF (myproc == 0) THEN
        WRITE (6, '(a)') 'Read in the Leaf Area Index data'
        READ (sfcunit,ERR=998) var2d
      END IF
      CALL mpisplit3d(var2d,nx,ny,1,lai)
    END IF

    IF(roufin == 1) THEN
      IF (myproc == 0) THEN
        WRITE (6, '(a)') 'Read in the surface roughness data'
        READ (sfcunit,ERR=998) var2d
      END IF
      CALL mpisplit3d(var2d,nx,ny,1,roufns)
    END IF

    IF(vegin == 1) THEN
      IF (myproc == 0) THEN
        WRITE (6, '(a)') 'Read in the vegetatin fraction data'
        READ (sfcunit,ERR=998) var2d
      END IF
      CALL mpisplit3d(var2d,nx,ny,1,veg)
    END IF

    IF (ndviin == 1) THEN
      IF (myproc == 0) THEN
        WRITE (6, '(a)') 'Read in the NDVI data'
        READ (sfcunit,ERR=998) ndvi
      END IF
      CALL mpisplit3d(var2d,nx,ny,1,ndvi)
    END IF

    IF (myproc == 0) THEN
      CLOSE ( sfcunit )
      CALL retunit ( sfcunit )
    END IF

  ELSE IF (sfcfmt == 3) THEN

    ALLOCATE(var2di(nxlg,nylg))
    ALLOCATE(var2d (nxlg,nylg))
    ALLOCATE(var3di(nxlg,nylg,nstyps))
    ALLOCATE(var3d (nxlg,nylg,nstyps))

    nstypin = MIN(nstyp1, nstyps)

    IF (myproc == 0) THEN
      CALL hdfrd3di(sd_id,"soiltyp",nxlg,nylg,nstypin,var3di,stat)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE (6, '(a)') 'Read in soiltyp'
        stypin = 1
      ELSE
        WRITE (6, '(a)') 'Variable soiltyp is not in the data set.'
        stypin = 0
      END IF
    END IF
    CALL mpisplit3di(var3di,nx,ny,nstypin,soiltyp)
    CALL mpupdatei(stypin,1)

    IF (myproc == 0) THEN
      CALL hdfrd3d(sd_id,"stypfrct",nxlg,nylg,nstypin,                &
                   var3d,stat,itmp,atmp1,atmp2)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE (6, '(a)') 'Read in stypfrct'
      ELSE
        WRITE (6, '(a)') 'Variable stypfrct is not in the data set.'
        var3d(:,:,1) = 1.
        var3d(:,:,2:nstypin) = 0.
      END IF
    END IF
    CALL mpisplit3d(var3d,nx,ny,nstypin,stypfrct)

    IF (myproc == 0) THEN
      CALL hdfrd2di(sd_id,"vegtyp",nxlg,nylg,var2di,stat)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE (6, '(a)') 'Read in vegtyp'
        vtypin = 1
      ELSE
        WRITE (6, '(a)') 'Variable vegtyp is not in the data set.'
        vtypin = 0
      END IF
    END IF
    CALL mpisplit2di(var2di,nx,ny,vegtyp)
    CALL mpupdatei(vtypin, 1)

    IF (myproc == 0) THEN
      CALL hdfrd2d(sd_id,"lai",nxlg,nylg,var2d,stat,itmp)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE (6, '(a)') 'Read in lai'
        laiin = 1
      ELSE
        WRITE (6, '(a)') 'Variable lai is not in the data set.'
        laiin = 0
      END IF
    END IF
    CALL mpisplit2d(var2d,nx,ny,lai)
    CALL mpupdatei(laiin, 1)

    IF (myproc == 0) THEN
      CALL hdfrd2d(sd_id,"roufns",nxlg,nylg,var2d,stat,itmp)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE (6, '(a)') 'Read in roufns'
        roufin = 1
      ELSE
        WRITE (6, '(a)') 'Variable roufns is not in the data set.'
        roufin = 0
      END IF
    END IF
    CALL mpisplit2d(var2d,nx,ny,roufns)
    CALL mpupdatei(roufin, 1)

    IF (myproc == 0) THEN
      CALL hdfrd2d(sd_id,"veg",nxlg,nylg,var2d,stat,itmp)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE (6, '(a)') 'Read in veg'
        vegin = 1
      ELSE
        WRITE (6, '(a)') 'Variable veg is not in the data set.'
        vegin = 0
      END IF
    END IF
    CALL mpisplit2d(var2d,nx,ny,veg)
    CALL mpupdatei(vegin, 1)

    IF (myproc == 0) THEN
      CALL hdfrd2d(sd_id,'ndvi',nxlg,nylg,var2d,stat,itmp)
      IF (stat > 1) GO TO 998
      IF (stat == 0) THEN
        WRITE (6, '(a)') 'Read in ndvi'
        ndviin = 1
      ELSE
        WRITE (6, '(a)') 'Variable ndvi is not in the data set.'
        ndviin = 0
      END IF
    END IF
    CALL mpisplit2d(var2d,nx,ny,ndvi)
    CALL mpupdatei(ndviin, 1)

    IF (myproc == 0) CALL hdfclose(sd_id, stat)

  ELSE IF (sfcfmt == 7) THEN    ! NetCDF format

    ALLOCATE(var2di(nxlg,nylg))
    ALLOCATE(var2d (nxlg,nylg))
    ALLOCATE(var3di(nxlg,nylg,nstyp1))
    ALLOCATE(var3d (nxlg,nylg,nstyp1))

    ALLOCATE(temi(nxlg-1,nylg-1,nstyp1))
    ALLOCATE(temr(nxlg-1,nylg-1,nstyp1))

    nstypin = MIN(nstyps,nstyp1)

    IF (myproc == 0) THEN

    WRITE (6, '(a/a,i2/a,i2/a,i2/a,i2/a,i2/a,i2)')                      &
                            ' Surface data flags for: ',                &
                            '        soiltyp --', stypin,               &
                            '         vegtyp --', vtypin,               &
                            '            lai --', laiin,                &
                            '         roufns --', roufin,               &
                            '            veg --', vegin,                &
                            '           ndvi --', ndviin

    WRITE (6, '(a/a,i2)') ' Number of soil types in each grid box:',    &
                          '          nstyp --', nstyp
    END IF


    IF(stypin == 1) THEN
      IF(myproc == 0) THEN
        CALL netread3di(sd_id,0,0,'SOILTYP', nxlg-1,nylg-1,nstyp1,temi)
        CALL netread3d (sd_id,0,0,'STYPFRCT',nxlg-1,nylg-1,nstyp1,temr)
        DO is = 1,nstypin
          DO j = 1,nylg-1
            DO i = 1,nxlg-1
              var3di(i,j,is) = temi(i,j,is)
              var3d (i,j,is) = temr(i,j,is)
            END DO
          END DO
        END DO
        CALL iedgfill(var3di,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nstypin,1,nstypin)
        CALL  edgfill(var3d, 1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,nstypin,1,nstypin)
      END IF
      CALL mpisplit3di(var3di,nx,ny,nstypin,soiltyp)
      CALL mpisplit3d (var3d,nx,ny,nstypin,stypfrct)
    END IF

    IF(vtypin == 1) THEN
      IF(myproc == 0) THEN
        CALL netread2di(sd_id,0,0,'VEGTYP',nxlg-1,nylg-1,temi)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var2di(i,j) = temi(i,j,1)
          END DO
        END DO
        CALL iedgfill(var2di,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2di(var2di,nx,ny,vegtyp)
    END IF

    IF(laiin == 1) THEN
      IF(myproc == 0) THEN
        CALL netread2d(sd_id,0,0,'LAI',nxlg-1,nylg-1,temr)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var2d(i,j) = temr(i,j,1)
          END DO
        END DO
        CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var2d,nx,ny,lai)
    END IF

    IF(roufin == 1) THEN
      IF(myproc == 0) THEN
        CALL netread2d(sd_id,0,0,'ROUFNS',nxlg-1,nylg-1,temr)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var2d(i,j) = temr(i,j,1)
          END DO
        END DO
        CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var2d,nx,ny,roufns)
    END IF

    IF(vegin == 1) THEN
      IF(myproc == 0) THEN
        CALL netread2d(sd_id,0,0,'VEG',nxlg-1,nylg-1,temr)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var2d(i,j) = temr(i,j,1)
          END DO
        END DO
        CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var2d,nx,ny,veg)
    END IF

    IF (ndviin == 1) THEN
      IF(myproc == 0) THEN
        CALL netread2d(sd_id,0,0,'NDVI',nxlg-1,nylg-1,var2d)
        DO j = 1,nylg-1
          DO i = 1,nxlg-1
            var2d(i,j) = temr(i,j,1)
          END DO
        END DO
        CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
      END IF
      CALL mpisplit2d(var2d,nx,ny,ndvi)
    END IF

    IF(myproc == 0) CALL netclose(sd_id)

    DEALLOCATE(temr, temi)

  END IF

  DEALLOCATE(var2d, var2di)
  DEALLOCATE(var3d, var3di)

  CALL fix_stypfrct_nstyp(nx,ny,nstyp1,nstyp,stypfrct)

  IF(stypin == 0) THEN
    WRITE(6,'(/1x,a,a)')                                                &
        'Soil type data not present in the surface propery ',           &
        'data file, set it to styp.'
    DO i=1,nx
      DO j=1,ny
        soiltyp(i,j,1) = styp
      END DO
    END DO
  END IF

  IF(vtypin == 0) THEN
    WRITE(6,'(/1x,a,a)')                                                &
        'Vegetation type data not present in the surface propery ',     &
        'data file, set it to vtyp.'

    DO i=1,nx-1
      DO j=1,ny-1
        vegtyp(i,j) = vtyp
      END DO
    END DO
  END IF

  IF(laiin == 0) THEN
    WRITE(6,'(/1x,a,a)')                                                &
        'Leaf Area Index data not present in the surface propery ',     &
        'data file, set it to lai0.'

    DO i=1,nx-1
      DO j=1,ny-1
        lai(i,j) = lai0
      END DO
    END DO
  END IF

  IF(roufin == 0) THEN
    WRITE(6,'(/1x,a,a)')                                                &
        'Roughness data not present in the surface propery ',           &
        'data file, set it to roufns0.'

    DO i=1,nx-1
      DO j=1,ny-1
        roufns(i,j) = roufns0
      END DO
    END DO
  END IF

  IF(vegin == 0) THEN
    WRITE(6,'(/1x,a,a)')                                                &
        'Vegetation fraction not present in the surface propery ',      &
        'data file, set it to veg0.'

    DO i=1,nx-1
      DO j=1,ny-1
        veg(i,j) = veg0
      END DO
    END DO
  END IF

  RETURN

  998   WRITE (6,'(/a,i2/a)')                                           &
         'READSPLITSFCDT: Read error in surface data file '                  &
         //sfcfile//' with the I/O unit ',sfcunit,                      &
         'The model will STOP in subroutine READSPLITSFCDT.'

  CALL arpsstop("arpsstop called from READSPLITSFCDT reading sfc file",1)

END SUBROUTINE readsplitsfcdt

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRITTRN                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE writtrn(nx,ny,ternfn,dx,dy,                                  &
           mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,         &
           hterain)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write out a terrain data file to be used by ARPS.
!
!  NOTE:
!
!    Any changes made here should also be made in subroutine writjointrn.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  5/1/1995.
!
!  MODIFICATION HISTORY:
!
!  2000/03/29 (Gene Bassett)
!  Removed the globcst.inc include.
!
!  2000/03/29 (Gene Bassett)
!  Added HDF4 format.
!
!  2004/08/10 (Yunheng Wang)
!  Added NetCDF 3.0 format.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!
!    mapproj  Type of map projection used to setup the analysis grid.
!    trulat1  1st real true latitude of map projection.
!    trulat2  2nd real true latitude of map projection.
!    trulon   Real true longitude of map projection.
!    sclfct   Map scale factor. At latitude = trulat1 and trulat2
!
!    dx       Model grid spacing in the x-direction east-west (meters)
!    dy       Model grid spacing in the y-direction east-west (meters)
!    ctrlat   Lat. at the origin of the model grid (deg. N)
!    ctrlon   Lon. at the origin of the model grid (deg. E)
!
!    ternfn   Terrain data file name
!    hterain  Terrain height (m)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER           :: nx,ny
  CHARACTER (LEN=*) :: ternfn

  REAL    :: dx, dy
  INTEGER :: mapproj       ! Map projection
  REAL    :: trulat1       ! 1st real true latitude of map projection
  REAL    :: trulat2       ! 2nd real true latitude of map projection
  REAL    :: trulon        ! Real true longitude of map projection
  REAL    :: sclfct        ! Map scale factor
  REAL    :: ctrlat        ! Center latitude of the model domain (deg. N)
  REAL    :: ctrlon        ! Center longitude of the model domain (deg. E)

  REAL    :: hterain(nx,ny)    ! The height of terrain.

  REAL, ALLOCATABLE :: var2d(:,:)
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: i,j
  INTEGER :: idummy, ierr, nunit
  REAL    :: rdummy

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
         ! unused array in hdf routines since NO COMPRESSION

  INTEGER :: istat, sd_id

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (terndmp == 0) RETURN

  WRITE (6,'(1x,/,2a)') 'WRITTRN: Opening Terrain file ',ternfn

!-----------------------------------------------------------------------
!
!  Write out in Fortran unformatted.
!
!-----------------------------------------------------------------------

  IF (terndmp == 1) THEN

    CALL getunit( nunit )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(ternfn, '-F f77 -N ieee', ierr)

    OPEN(nunit,FILE=trim(ternfn),FORM='unformatted',STATUS='unknown')

    WRITE(nunit) nx,ny

    idummy = 0
    rdummy = 0.0

    WRITE(nunit) idummy,mapproj,idummy,idummy,idummy,                &
                 idummy,idummy,idummy,idummy,idummy,                 &
                 idummy,idummy,idummy,idummy,idummy,                 &
                 idummy,idummy,idummy,idummy,idummy

    WRITE(nunit) dx   ,dy   ,ctrlat,ctrlon,rdummy ,                  &
                 rdummy,trulat1,trulat2,trulon,sclfct,               &
                 rdummy,rdummy,rdummy,rdummy,rdummy,                 &
                 rdummy,rdummy,rdummy,rdummy,rdummy

    WRITE(nunit) hterain

    CLOSE(nunit)

  ELSE IF (terndmp == 3) THEN

!-----------------------------------------------------------------------
!
!  Write out in HDF4.
!
!-----------------------------------------------------------------------

    CALL hdfopen(trim(ternfn), 2, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,'(a,a)') "WRITTRN: ERROR creating HDF4 file: ",          &
                  trim(ternfn)
      CALL arpsstop('ERROR: creating terrain HDF4 file.',1)
    END IF

    CALL hdfwrti(sd_id, 'nx', nx, istat)
    CALL hdfwrti(sd_id, 'ny', ny, istat)
    CALL hdfwrtr(sd_id, 'dx', dx, istat)
    CALL hdfwrtr(sd_id, 'dy', dy, istat)
    CALL hdfwrti(sd_id, 'mapproj', mapproj, istat)
    CALL hdfwrtr(sd_id, 'trulat1', trulat1, istat)
    CALL hdfwrtr(sd_id, 'trulat2', trulat2, istat)
    CALL hdfwrtr(sd_id, 'trulon',  trulon,  istat)
    CALL hdfwrtr(sd_id, 'sclfct',  sclfct,  istat)
    CALL hdfwrtr(sd_id, 'ctrlat',  ctrlat,  istat)
    CALL hdfwrtr(sd_id, 'ctrlon',  ctrlon,  istat)

    CALL hdfwrt2d(hterain,nx,ny,sd_id,0,0,                              &
                  'hterain','Terrain Height','m',itmp)
    CALL hdfclose(sd_id,istat)
    IF (istat /= 0) THEN
      WRITE (6,'(a,a)') "WRITTRN: ERROR on closing file ",trim(ternfn), &
                  " (status",istat,")"
    END IF
    !wdt end block

!-----------------------------------------------------------------------
!
!  Write out in NetCDF 3.0 format
!
!-----------------------------------------------------------------------
  ELSE IF (terndmp == 7) THEN

    ALLOCATE(var2d(nx-1,ny-1), STAT = istat)

!-----------------------------------------------------------------------
!
!  Define terrain file dimension and variables
!
!-----------------------------------------------------------------------

    CALL netopen(TRIM(ternfn), 'C', sd_id)

    CALL net_define_trn(sd_id,nx,ny,dx,dy,mapproj,sclfct,             &
                        trulat1,trulat2,trulon,ctrlat,ctrlon,istat)

    DO j = 1,ny-1
      DO i = 1, nx-1
        var2d(i,j) = hterain(i,j)
      END DO
    END DO

    CALL netwrt2d(sd_id,0,0,'HTERAIN',var2d,nx-1,ny-1)

    CALL netclose(sd_id)

    DEALLOCATE(var2d)

  ELSE

    ! alternate dump format ...
    WRITE(6,*) 'The supported terrain data dump format are ',           &
               'binary (terndmp=1), HDF4 no compressed (terndmp = 3),', &
               'and NetCDF 3 (terndmp = 7).'
    CALL arpsstop('Terrain data dump format is not supported.',1)

  END IF

  RETURN
END SUBROUTINE writtrn
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRITJOINTRN                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE writjointrn(nx,ny,ternfn,dx,dy,                              &
           mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,         &
           hterain)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Join and Write out a terrain data file to be used by ARPS.
!
!  NOTE:
!  Any changes made here should also be made in subroutine writtrn.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/20/2003
!  Modified from writtrn.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    mapproj  Type of map projection used to setup the analysis grid.
!    trulat1  1st real true latitude of map projection.
!    trulat2  2nd real true latitude of map projection.
!    trulon   Real true longitude of map projection.
!    sclfct   Map scale factor. At latitude = trulat1 and trulat2
!
!    dx       Model grid spacing in the x-direction east-west (meters)
!    dy       Model grid spacing in the y-direction east-west (meters)
!    ctrlat   Lat. at the origin of the model grid (deg. N)
!    ctrlon   Lon. at the origin of the model grid (deg. E)
!
!    ternfn   Terrain data file name
!    hterain  Terrain height (m)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER           :: nx,ny
  CHARACTER (LEN=*) :: ternfn

  REAL    :: dx
  REAL    :: dy
  INTEGER :: mapproj          ! Map projection
  REAL    :: trulat1          ! 1st real true latitude of map projection
  REAL    :: trulat2          ! 2nd real true latitude of map projection
  REAL    :: trulon           ! Real true longitude of map projection
  REAL    :: sclfct           ! Map scale factor
  REAL    :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL    :: ctrlon           ! Center longitude of the model domain (deg. E)

  REAL    :: hterain(nx,ny)   ! The height of terrain.

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: i,j
  INTEGER :: idummy, ierr, nunit
  REAL    :: rdummy

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
         ! unused array in hdf routines since NO COMPRESSION

  INTEGER :: istat, sd_id

  INTEGER           :: nxlg, nylg
  REAL, ALLOCATABLE :: ht_global(:,:)
  REAL, ALLOCATABLE :: var2d(:,:)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (terndmp == 0) RETURN

  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

  ALLOCATE(ht_global(nxlg,nylg), STAT = istat)

  IF ( myproc == 0 )   &
    WRITE (6,'(a,a)') 'WRITJOINTRN: Opening Terrain file ',ternfn

  CALL mpimerge2d(hterain,nx,ny,ht_global)

  IF (myproc ==0 ) THEN

!-----------------------------------------------------------------------
!
!  Write out in Fortran unformatted.
!
!-----------------------------------------------------------------------

    IF (terndmp == 1) THEN

      CALL getunit( nunit )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(ternfn, '-F f77 -N ieee', ierr)

      OPEN(nunit,FILE=trim(ternfn),FORM='unformatted',STATUS='unknown')

      WRITE(nunit) nxlg,nylg

      idummy = 0
      rdummy = 0.0

      WRITE(nunit) idummy,mapproj,idummy,idummy,idummy,                &
                   idummy,idummy,idummy,idummy,idummy,                 &
                   idummy,idummy,idummy,idummy,idummy,                 &
                   idummy,idummy,idummy,idummy,idummy

      WRITE(nunit) dx   ,dy   ,ctrlat,ctrlon,rdummy ,                  &
                   rdummy,trulat1,trulat2,trulon,sclfct,               &
                   rdummy,rdummy,rdummy,rdummy,rdummy,                 &
                   rdummy,rdummy,rdummy,rdummy,rdummy

      WRITE(nunit) ht_global      ! only root process write

      CLOSE(nunit)

!-----------------------------------------------------------------------
!
!  Write out in HDF4.
!
!-----------------------------------------------------------------------

    ELSE IF (terndmp == 3) THEN

      CALL hdfopen(trim(ternfn), 2, sd_id)
      IF (sd_id < 0) THEN
        WRITE (6,'(a,a)') "WRITJOINTRN: ERROR creating HDF4 file: ",    &
                          trim(ternfn)
        CALL arpsstop('ERROR: creating terrain HDF4 file.',1)
      END IF

      CALL hdfwrti(sd_id, 'nx', nxlg, istat)
      CALL hdfwrti(sd_id, 'ny', nylg, istat)
      CALL hdfwrtr(sd_id, 'dx', dx,   istat)
      CALL hdfwrtr(sd_id, 'dy', dy,   istat)
      CALL hdfwrti(sd_id, 'mapproj', mapproj, istat)
      CALL hdfwrtr(sd_id, 'trulat1', trulat1, istat)
      CALL hdfwrtr(sd_id, 'trulat2', trulat2, istat)
      CALL hdfwrtr(sd_id, 'trulon',  trulon,  istat)
      CALL hdfwrtr(sd_id, 'sclfct',  sclfct,  istat)
      CALL hdfwrtr(sd_id, 'ctrlat',  ctrlat,  istat)
      CALL hdfwrtr(sd_id, 'ctrlon',  ctrlon,  istat)

      CALL hdfwrt2d(ht_global,nxlg,nylg,sd_id,0,0,                      &
                    'hterain','Terrain Height','m',itmp)

      CALL hdfclose(sd_id,istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "WRITJOINTRN: ERROR on closing file ",trim(ternfn), &
                    " (status",istat,")"
      END IF

!-----------------------------------------------------------------------
!
!  Write out in NetCDF 3.0 format
!
!-----------------------------------------------------------------------

    ELSE IF (terndmp == 7) THEN

      ALLOCATE(var2d(nxlg-1,nylg-1), STAT = istat)

!-----------------------------------------------------------------------
!
!  Define terrain file dimension and variables
!
!-----------------------------------------------------------------------

      CALL netopen(TRIM(ternfn), 'C', sd_id)
      CALL net_define_trn(sd_id,nxlg,nylg,dx,dy,mapproj,sclfct,         &
                          trulat1,trulat2,trulon,ctrlat,ctrlon,istat)

      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var2d(i,j) = ht_global(i,j)
        END DO
      END DO
      CALL netwrt2d(sd_id,0,0,'HTERAIN',var2d,nxlg-1,nylg-1)

      CALL netclose(sd_id)

      DEALLOCATE(var2d)

    ELSE

      ! alternate dump format ...
      WRITE(6,*) 'The supported terrain data dump format are ',         &
               'binary (terndmp=1), HDF4 no compressed (terndmp = 3),', &
               'and NetCDF 3 (terndmp = 7).'
      CALL arpsstop('Terrain data dump format is not supported.',1)

    END IF      ! terndmp

  END IF  ! myproc == 0

  DEALLOCATE(ht_global)

  RETURN
END SUBROUTINE writjointrn
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READTRN                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readtrn(nx,ny, dx,dy, ternfile,                              &
                   mapproj,trulat1,trulat2,trulon,sclfct,               &
                   ctrlat,ctrlon,hterain )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the terrain data into model array hterain from a specified
!  terrain data file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  2/27/1994.
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  2000/03/29 (Gene Bassett)
!  Removed the globcst.inc include.
!
!  2000/03/29 (Gene Bassett)
!  Added HDF4 format.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    dx       Grid interval in x-direction
!    dy       Grid interval in y-direction
!
!    ternfile    Terrain data file name
!
!  mapproj    Type of map projection used to setup the analysis grid.
!  trulat1    1st real true latitude of map projection.
!  trulat2    2nd real true latitude of map projection.
!  trulon     Real true longitude of map projection.
!  sclfct     Map scale factor. At latitude = trulat1 and trulat2
!
!  ctrlat    Lat. at the origin of the model grid (deg. N)
!  ctrlon    Lon. at the origin of the model grid (deg. E)
!
!  OUTPUT:
!
!    hterain  Terrain height (m)
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

  INTEGER :: nx            ! Number of grid points in the x-direction
  INTEGER :: ny            ! Number of grid points in the y-direction
  REAL    :: dx            ! Grid interval in x-direction
  REAL    :: dy            ! Grid interval in y-direction
  CHARACTER(LEN=*) :: ternfile ! Terrain data file name

  INTEGER :: mapproj       ! Map projection
  REAL :: trulat1          ! 1st real true latitude of map projection
  REAL :: trulat2          ! 2nd real true latitude of map projection
  REAL :: trulon           ! Real true longitude of map projection
  REAL :: sclfct           ! Map scale factor
  REAL :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon           ! Center longitude of the model domain (deg. E)

  REAL :: hterain(nx,ny)   ! Terrain height.

  REAL, ALLOCATABLE :: var2d(:,:)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  INTEGER :: inunit,istat
  INTEGER :: nxin,nyin,idummy,ierr
  REAL    :: dxin,dyin,rdummy,amin,amax

  INTEGER :: ireturn
  REAL    :: ctrlatin,ctrlonin,                                         &
             trulat1in,trulat2in,trulonin,sclfctin
  INTEGER :: mapprojin

  CHARACTER (LEN=256) :: savename

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
         ! unused array in hdf routines since NO COMPRESSION

  INTEGER :: sd_id

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Read in the terrain data.
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    savename(1:256) = ternfile(1:256)
    CALL gtsplitfn(savename,1,1,loc_x,loc_y,1,1,0,0,1,lvldbg,ternfile,ireturn)
  END IF

  WRITE (6,'(1x,2a)') 'READTRN: reading in external supplied terrain '// &
                      'data from file ',trim(ternfile)

!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------

  IF (ternfmt == 1) THEN

    CALL getunit( inunit )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(ternfile, '-F f77 -N ieee', ierr)

    OPEN(UNIT=inunit,FILE=trim(ternfile),FORM='unformatted',            &
         STATUS='old',IOSTAT=istat)

    IF( istat /= 0) THEN

      WRITE(6,'(/1x,a,a,/1x,a/)')                                       &
          'Error occured when opening terrain data file ',              &
          ternfile,' Job stopped in READTRN.'
      CALL arpsstop("arpsstop called from READTRN reading terrain file",1)

    END IF

    READ(inunit,ERR=999) nxin,nyin

    READ(inunit,ERR=999) idummy,idummy,idummy,idummy,idummy,            &
               idummy,idummy,idummy,idummy,idummy,                      &
               idummy,idummy,idummy,idummy,idummy,                      &
               idummy,idummy,idummy,idummy,idummy

    READ(inunit,ERR=999) dxin  ,dyin  ,rdummy,rdummy,rdummy,            &
               rdummy,rdummy,rdummy,rdummy,rdummy,                      &
               rdummy,rdummy,rdummy,rdummy,rdummy,                      &
               rdummy,rdummy,rdummy,rdummy,rdummy

  ELSE IF(ternfmt == 3) THEN

    CALL hdfopen(trim(ternfile), 1, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "READTRN: ERROR opening ",                            &
                 trim(ternfile)," as HDF 4 format for reading."
      GO TO 999
    END IF

    CALL hdfrdi(sd_id,"nx",nxin,istat)
    CALL hdfrdi(sd_id,"ny",nyin,istat)
    CALL hdfrdr(sd_id,"dx",dxin,istat)
    CALL hdfrdr(sd_id,"dy",dyin,istat)
    CALL hdfrdi(sd_id,"mapproj",mapprojin,istat)
    CALL hdfrdr(sd_id,"trulat1",trulat1in,istat)
    CALL hdfrdr(sd_id,"trulat2",trulat2in,istat)
    CALL hdfrdr(sd_id,"trulon",trulonin,istat)
    CALL hdfrdr(sd_id,"sclfct",sclfctin,istat)
    CALL hdfrdr(sd_id,"ctrlat",ctrlatin,istat)
    CALL hdfrdr(sd_id,"ctrlon",ctrlonin,istat)
!wdt end block

  ELSE IF (ternfmt == 7) THEN

    CALL netopen(TRIM(ternfile),'R',sd_id)

    CALL net_get_trn(sd_id,nxin,nyin,dxin,dyin,mapprojin,sclfctin,      &
                  trulat1in,trulat2in,trulonin,ctrlatin,ctrlonin, istat)

  ELSE

    ! alternate dump format ...
    WRITE(6,*) 'The supported terrain data format are ',                &
               'binary (ternfmt=1) and HDF4 no compressed (ternfmt = 3).'
    CALL arpsstop('Terrain data format is not supported.',1)

  END IF

  IF (ternfmt == 1) THEN
    WRITE (6,*)                                                         &
        "READTRN: WARNING, not checking map projection parameters"
    CALL checkgrid2d(nx,ny,nxin,nyin,                                   &
              dx,dy,ctrlat,ctrlon,                                      &
              mapproj,trulat1,trulat2,trulon,sclfct,                    &
              dxin,dyin,ctrlat,ctrlon,                                  &
              mapproj,trulat1,trulat2,trulon,sclfct,ireturn)
  ELSE
    CALL checkgrid2d(nx,ny,nxin,nyin,                                   &
             dx,dy,ctrlat,ctrlon,                                       &
             mapproj,trulat1,trulat2,trulon,sclfct,                     &
             dxin,dyin,ctrlatin,ctrlonin,                               &
             mapprojin,trulat1in,trulat2in,trulonin,sclfctin,ireturn)
  END IF

  IF (ireturn /= 0) THEN
    WRITE (6,*) "READTRN: ERROR, grid parameter mismatch"
    CALL arpsstop("arpsstop called from READTRN parameter mismatch",1)
  END IF

  IF (ternfmt == 1) THEN

    READ(inunit,ERR=999) hterain

    CALL retunit( inunit )
    CLOSE (UNIT=inunit)

  ELSE IF (ternfmt == 3) THEN

    CALL hdfrd2d(sd_id,"hterain",nx,ny,hterain,istat,itmp)
    CALL hdfclose(sd_id,istat)

  ELSE IF (ternfmt == 7) THEN

    ALLOCATE(var2d(nx-1,ny-1), STAT = istat)

    CALL netread2d(sd_id,0,0,'HTERAIN',nx-1,ny-1,var2d)
    CALL netclose(sd_id)
    DO j = 1,ny-1
      DO i = 1,nx-1
        hterain(i,j) = var2d(i,j)
      END DO
    END DO
    CALL edgfill(hterain,1,nx, 1, nx-1, 1,ny, 1,ny-1, 1,1,1,1)

    DEALLOCATE(var2d)

  END IF

  WRITE(6,'(3x,a)') 'Minimum and maximum terrain height:'

  CALL a3dmax0(hterain,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1,amax,amin)
  WRITE(6,'(3x,2(a,e13.6),/)') 'htrnmin = ', amin,', htrnmax=',amax

  IF (mp_opt > 0) ternfile(1:256) = savename(1:256)

  RETURN

  999   WRITE(6,'(1x,a)')                                               &
        'Error in reading terrain data. Job stopped in READTRN.'
  CALL arpsstop("arpsstop called from READTRN reading file",1)

END SUBROUTINE readtrn
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READSPLITTRN               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readsplittrn(nx,ny, dx,dy, ternfile,                         &
           mapproj,trulat1,trulat2,trulon,sclfct,                       &
           ctrlat,ctrlon,hterain )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the terrain data into model array hterain from a specified
!  terrain data file, and split and scatter to each processors for
!  MPI processes.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  8/30/2002.
!  Based on subroutine readtrn
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    dx       Grid interval in x-direction
!    dy       Grid interval in y-direction
!
!    ternfile    Terrain data file name
!
!  mapproj    Type of map projection used to setup the analysis grid.
!  trulat1    1st real true latitude of map projection.
!  trulat2    2nd real true latitude of map projection.
!  trulon     Real true longitude of map projection.
!  sclfct     Map scale factor. At latitude = trulat1 and trulat2
!
!  ctrlat    Lat. at the origin of the model grid (deg. N)
!  ctrlon    Lon. at the origin of the model grid (deg. E)
!
!  OUTPUT:
!
!    hterain  Terrain height (m)
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

  INTEGER :: nx            ! Number of grid points in the x-direction
  INTEGER :: ny            ! Number of grid points in the y-direction
  REAL    :: dx            ! Grid interval in x-direction
  REAL    :: dy            ! Grid interval in y-direction
  CHARACTER (LEN=*) :: ternfile ! Terrain data file name

  INTEGER :: mapproj       ! Map projection
  REAL :: trulat1          ! 1st real true latitude of map projection
  REAL :: trulat2          ! 2nd real true latitude of map projection
  REAL :: trulon           ! Real true longitude of map projection
  REAL :: sclfct           ! Map scale factor
  REAL :: ctrlat           ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon           ! Center longitude of the model domain (deg. E)

  REAL :: hterain(nx,ny)   ! Terrain height.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  INTEGER :: inunit,istat
  INTEGER :: nxin,nyin,idummy,ierr
  REAL :: dxin,dyin,rdummy,amin,amax

  INTEGER :: ireturn
  REAL :: ctrlatin,ctrlonin,                                            &
          trulat1in,trulat2in,trulonin,sclfctin
  INTEGER :: mapprojin

  INTEGER(KIND=selected_int_kind(4)) :: itmp(1)
         ! unused array in hdf routines since NO COMPRESSION

  INTEGER :: sd_id

  INTEGER :: nxlg, nylg, old_mp_opt
  REAL, ALLOCATABLE :: tem(:,:)
  REAL, ALLOCATABLE :: var2d(:,:)

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3

  ALLOCATE(tem(nxlg, nylg), stat= istat)
!
!-----------------------------------------------------------------------
!
!  Read in header information.
!
!-----------------------------------------------------------------------

  IF (myproc == 0) THEN

    WRITE (6,*) "READSPLITTRN: reading in external supplied terrain ",  &
                "data from file ",trim(ternfile)

    IF (ternfmt == 1) THEN

      CALL getunit( inunit )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(ternfile, '-F f77 -N ieee', ierr)

      OPEN(UNIT=inunit,FILE=trim(ternfile),FORM='unformatted',            &
           STATUS='old',IOSTAT=istat)

      IF( istat /= 0) THEN

        WRITE(6,'(/1x,a,a,/1x,a/)')                                       &
            'Error occured when opening terrain data file ',              &
            ternfile,' Job stopped in READSPLITTRN.'
        CALL arpsstop("arpsstop called from READSPLITTRN reading terrain file",1)

      END IF

      READ(inunit,ERR=999) nxin,nyin

      READ(inunit,ERR=999) idummy,idummy,idummy,idummy,idummy,            &
                 idummy,idummy,idummy,idummy,idummy,                      &
                 idummy,idummy,idummy,idummy,idummy,                      &
                 idummy,idummy,idummy,idummy,idummy

      READ(inunit,ERR=999) dxin  ,dyin  ,rdummy,rdummy,rdummy,            &
                 rdummy,rdummy,rdummy,rdummy,rdummy,                      &
                 rdummy,rdummy,rdummy,rdummy,rdummy,                      &
                 rdummy,rdummy,rdummy,rdummy,rdummy

    ELSE IF (ternfmt == 3) THEN

      CALL hdfopen(trim(ternfile), 1, sd_id)
      IF (sd_id < 0) THEN
        WRITE (6,*) "READSPLITTRN: ERROR opening ",                       &
                   trim(ternfile)," for reading."
        GO TO 999
      END IF

      CALL hdfrdi(sd_id,"nx",nxin,istat)
      CALL hdfrdi(sd_id,"ny",nyin,istat)
      CALL hdfrdr(sd_id,"dx",dxin,istat)
      CALL hdfrdr(sd_id,"dy",dyin,istat)
      CALL hdfrdi(sd_id,"mapproj",mapprojin,istat)
      CALL hdfrdr(sd_id,"trulat1",trulat1in,istat)
      CALL hdfrdr(sd_id,"trulat2",trulat2in,istat)
      CALL hdfrdr(sd_id,"trulon",trulonin,istat)
      CALL hdfrdr(sd_id,"sclfct",sclfctin,istat)
      CALL hdfrdr(sd_id,"ctrlat",ctrlatin,istat)
      CALL hdfrdr(sd_id,"ctrlon",ctrlonin,istat)

    ELSE IF (ternfmt == 7) THEN

      CALL netopen(TRIM(ternfile), 'R', sd_id)
      CALL net_get_trn(sd_id,nxin,nyin,dxin,dyin,mapprojin,sclfctin,    &
                  trulat1in,trulat2in,trulonin,ctrlatin,ctrlonin, istat)

    ELSE

      ! alternate dump format ...
      WRITE(6,*) 'The supported terrain data format are ',              &
             'binary (ternfmt=1) and HDF4 no compressed (ternfmt = 3).'
      CALL arpsstop('Terrain data format is not supported.',1)

    END IF

    nxin = (nxin - 3)/nproc_x + 3
    nyin = (nyin - 3)/nproc_y + 3

  END IF

  CALL mpupdatei(nxin,1)
  CALL mpupdatei(nyin,1)
  CALL mpupdater(dxin,1)
  CALL mpupdater(dyin,1)

  IF (ternfmt == 1) THEN
    WRITE (6,*)                                                         &
        "READSPLITTRN: WARNING, not checking map projection parameters"
    CALL checkgrid2d(nx,ny,nxin,nyin,                                   &
              dx,dy,ctrlat,ctrlon,                                      &
              mapproj,trulat1,trulat2,trulon,sclfct,                    &
              dxin,dyin,ctrlat,ctrlon,                                  &
              mapproj,trulat1,trulat2,trulon,sclfct,ireturn)
  ELSE
    CALL mpupdatei(mapprojin,1)
    CALL mpupdater(trulat1in,1)
    CALL mpupdater(trulat2in,1)
    CALL mpupdater(trulonin, 1)
    CALL mpupdater(sclfctin, 1)
    CALL mpupdater(ctrlatin, 1)
    CALL mpupdater(ctrlonin, 1)

    CALL checkgrid2d(nx,ny,nxin,nyin,                                   &
             dx,dy,ctrlat,ctrlon,                                       &
             mapproj,trulat1,trulat2,trulon,sclfct,                     &
             dxin,dyin,ctrlatin,ctrlonin,                               &
             mapprojin,trulat1in,trulat2in,trulonin,sclfctin,ireturn)
  END IF

  IF (ireturn /= 0) THEN
    WRITE (6,*) "READSPLITTRN: ERROR, grid parameter mismatch"
    CALL arpsstop("arpsstop called from READSPLITTRN parameter mismatch",1)
  END IF

  IF (myproc == 0) THEN

    IF (ternfmt == 1) THEN

      READ(inunit,ERR=999) tem

      CALL retunit( inunit )
      CLOSE (UNIT=inunit)

    ELSE IF (ternfmt == 3) THEN

      CALL hdfrd2d(sd_id,"hterain",nxlg,nylg,tem,istat,itmp)

      CALL hdfclose(sd_id,istat)

    ELSE IF (ternfmt == 7) THEN

      ALLOCATE(var2d(nxlg-1,nylg-1), STAT = istat)
      CALL netread2d(sd_id,0,0,'HTERAIN',nxlg-1,nylg-1,var2d)
      CALL netclose(sd_id)
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          tem(i,j) = var2d(i,j)
        END DO
      END DO
      CALL edgfill(tem,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)

      DEALLOCATE(var2d)

    END IF

    old_mp_opt = mp_opt
    mp_opt = 0

    WRITE(6,'(1x,a/)') 'Minimum and maximum terrain height:'
    CALL a3dmax0(tem,1,nxlg,1,nxlg-1, 1,nylg,1,nylg-1, 1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'htrnmin = ', amin,', htrnmax=',amax

    mp_opt = old_mp_opt

  END IF

  CALL mpisplit2d(tem,nx,ny,hterain)

  DEALLOCATE(tem)

  RETURN

  999   WRITE(6,'(1x,a)')                                               &
        'Error in reading terrain data. Job stopped in READSPLITTRN.'
  CALL arpsstop("arpsstop called from READSPLITTRN reading file",1)

END SUBROUTINE readsplittrn

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRITESND                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE writesnd(nx,ny,nz,                                           &
           ubar,vbar,ptbar,pbar,qvbar,                                  &
           zp, rhobar, zs)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Print the base state sounding profile at the center of the domain.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  5/1/98 (G. Bassett, D. Weber)
!  Moved from inibase.
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
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    zp       Vertical coordinate of grid points in physical space(m)
!
!    rhobar   Base state air density (kg/m**3)
!    zs       Physical coordinate height at scalar point (m)
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

  INTEGER :: nx,ny,nz          ! The number of grid points in 3
                               ! directions

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: qvbar (nx,ny,nz)     ! Base state water specific humidity
                               ! (kg/kg).
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of staggered grid.

  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: zs    (nx,ny,nz)     ! Physical coordinate height at scalar
                               ! point (m)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k, istat, nunit, itrnmin, jtrnmin
  REAL :: ternmin
  CHARACTER(LEN=256) :: outdirname
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(myproc == 0) THEN

    ternmin = zp(1,1,2)
    itrnmin = 1
    jtrnmin = 1

    DO j=1,ny-1
      DO i=1,nx-1
        IF(zp(i,j,2) < ternmin) THEN
          ternmin = zp(i,j,2)
          itrnmin = i
          jtrnmin = j
        END IF
      END DO
    END DO

    i = itrnmin
    j = jtrnmin

    WRITE(6,'(/a,i4,a,i4/)')                                            &
        ' Base state z-profiles at I=',i,', J=',j

    CALL getunit( nunit )

    CALL get_output_dirname(0,dirname,curtim,1,outdirname,istatus)

    OPEN(UNIT=nunit, FILE=TRIM(outdirname)//runname(1:lfnkey)//'.sound',&
         STATUS='unknown',IOSTAT=istat)

    IF(istat /= 0) THEN
      WRITE (6,'(1x,a,a)') 'Error opening file ',                       &
          runname(1:lfnkey)//'.sound'
      CALL arpsstop("arpsstop called from WRITESND opening file",1)
    END IF

    WRITE(6,'(1x,2a)')                                                  &
        '   k      z(m)    pbar(Pascal)  ptbar(K) rhobar(kg/m**3)',     &
        ' Qvbar(kg/kg)  U(m/s)     V(m/s)'
    WRITE(nunit,'(1x,2a)')                                              &
        '   k      z(m)    pbar(Pascal)  ptbar(K) rhobar(kg/m**3)',     &
        ' Qvbar(kg/kg)  U(m/s)     V(m/s)'


    DO k=nz-1,1,-1

      WRITE(6, '(1x,i4,2f12.3,f12.5,2f12.8,2f12.5)')                    &
          k,zs(i,j,k),pbar(i,j,k),ptbar(i,j,k),rhobar(i,j,k)            &
          ,qvbar(i,j,k),ubar(i,j,k), vbar(i,j,k)

      WRITE(nunit,'(1x,i4,2f12.3,f12.5,2f12.8,2f12.5)')                 &
          k,zs(i,j,k),pbar(i,j,k),ptbar(i,j,k),rhobar(i,j,k)            &
          ,qvbar(i,j,k),ubar(i,j,k), vbar(i,j,k)

    END DO

    CLOSE( UNIT = nunit )
    CALL retunit( nunit )

  END IF     !Message Passing

  RETURN
END SUBROUTINE writesnd
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE SFCCNTL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sfccntl(nx,ny, sfcoutfl,                                     &
           stypout,vtypout,laiout,rfnsout,vegout,ndviout,               &
           x,y, temxy1,temxy2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  ARPS surface property data description file for GrADS display
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  10/05/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    x        The x-coord. of the physical and computational grid.
!             Defined at u-point.
!    y        The y-coord. of the physical and computational grid.
!             Defined at v-point.
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny             ! Number of grid points in 3 directions

  CHARACTER (LEN=*   ) :: sfcoutfl ! Surface data file name

  INTEGER :: stypout           ! Flag for output of soiltyp
  INTEGER :: vtypout           ! Flag for output of vegtyp
  INTEGER :: laiout            ! Flag for output of lai
  INTEGER :: rfnsout           ! Flag for output of roufns
  INTEGER :: vegout            ! Flag for output of veg
  INTEGER :: ndviout           ! Flag for output of ndvi

  REAL :: x(nx)                ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y(ny)                ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.

  REAL :: temxy1(nx,ny)        ! Temporary array
  REAL :: temxy2(nx,ny)        ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: varnumax
  PARAMETER ( varnumax = 100 )

  CHARACTER (LEN=256) :: sfcctlfl
  CHARACTER (LEN=256) :: sfcdatfl
  CHARACTER (LEN=15)  :: chrstr, chr1
  CHARACTER (LEN=8)   :: varnam(varnumax)
  CHARACTER (LEN=60)  :: vartit(varnumax)
  CHARACTER (LEN=10)  :: varparam(varnumax)

  INTEGER :: varlev(varnumax)

  INTEGER :: varnum
  INTEGER :: nchout0

  CHARACTER (LEN=3) :: monnam(12)            ! Name of months
  DATA monnam/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                 &
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/

  CHARACTER (LEN=2) :: dtunit
  INTEGER :: ntm,tinc

  REAL :: lonmin,latmin,lonmax,latmax
  REAL :: xbgn,ybgn
  REAL :: xinc,yinc
  REAL :: lat11,lat12,lon11,lon12,lat21,lat22,lon21,lon22
  REAL :: latinc,loninc

  INTEGER :: ctllen,fnlen
  INTEGER :: i,k, is
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
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
!  Open the GrADS data control file: runname.ctl
!
!-----------------------------------------------------------------------
!
  fnlen = LEN( sfcoutfl )
  CALL strlnth( sfcoutfl, fnlen )

  ctllen = fnlen + 4
  sfcctlfl(1:ctllen) = sfcoutfl(1:fnlen)//'.ctl'

  CALL fnversn( sfcctlfl, ctllen )
  CALL getunit (nchout0)
  WRITE (6,'(a)') 'The GrADS data control file is '                     &
                    //sfcctlfl(1:ctllen)

  OPEN (nchout0, FILE = sfcctlfl(1:ctllen), STATUS = 'unknown')

  xbgn = 0.5 * (x(1) + x(2))
  ybgn = 0.5 * (y(1) + y(2))

  xinc = (x(2) - x(1))
  yinc = (y(2) - y(1))

  CALL xytoll(nx,ny,x,y,temxy1,temxy2)

  CALL xytoll(1,1,xbgn,ybgn,lat11,lon11)

  CALL a3dmax0(temxy1,1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,                    &
               latmax,latmin)
  CALL a3dmax0(temxy2,1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,                    &
               lonmax,lonmin)

  latinc = (latmax-latmin)/(ny-1)
  loninc = (lonmax-lonmin)/(nx-1)

  WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                                 &
           'latmin:latmax:latinc = ',                                   &
            latmin,':',latmax,':',latinc
  WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                                 &
           'lonmin:lonmax:loninc = ',                                   &
           lonmin,':',lonmax,':',loninc

  WRITE (chrstr,'(i2.2,a,i2.2,a,i2.2,a3,i4.4)')                         &
      hour,':',minute,'Z',day,monnam(month),year

  ntm    = 1
  tinc   = 1
  dtunit = 'MO'

  DO i=1,varnumax
    varlev(i) = 0
    varparam(i) = '99'
  END DO

  varnum = 0

  IF ( stypout == 1 ) THEN
    IF ( nstyp == 1 ) THEN
      varnum = varnum + 1
      varnam(varnum) = 'styp'
      vartit(varnum) = 'Soil type'
      varparam(varnum) = '-1,40,4'
    ELSE
      DO is=1,nstyp
        WRITE (chr1,'(i2.2)') is

        varnum = varnum + 1
        varnam(varnum) = 'styp'//chr1(1:2)
        vartit(varnum) = 'Soil type '//chr1(1:2)
        varparam(varnum) = '-1,40,4'

        varnum = varnum + 1
        varnam(varnum) = 'sfct'//chr1(1:2)
        vartit(varnum) = 'Fraction of soil type '//chr1(1:2)
      END DO
    END IF
  END IF

  IF ( vtypout == 1 ) THEN
    varnum = varnum + 1
    varnam(varnum) = 'vtyp'
    vartit(varnum) = 'Vegetation type'
    varparam(varnum) = '-1,40,4'
  END IF

  IF ( laiout == 1 ) THEN
    varnum = varnum + 1
    varnam(varnum) = 'lai'
    vartit(varnum) = 'Leaf Area Index'
  END IF

  IF ( rfnsout == 1 ) THEN
    varnum = varnum + 1
    varnam(varnum) = 'roufns'
    vartit(varnum) = 'Surface roughness (m)'
  END IF

  IF ( vegout == 1 ) THEN
    varnum = varnum + 1
    varnam(varnum) = 'veg'
    vartit(varnum) = 'Vegetation fraction'
  END IF

  IF ( ndviout == 1 ) THEN
    varnum = varnum + 1
    varnam(varnum) = 'ndvi'
    vartit(varnum) = 'NDVI'
  END IF

  WRITE (nchout0,'(a/a)')                                               &
      'TITLE   ARPS surface property data control for '                 &
      //runname(1:lfnkey),'*'

  WRITE (nchout0,'(a,a)')                                               &
      'DSET    ',sfcoutfl(1:fnlen)

  WRITE (nchout0,'(a)')                                                 &
      'OPTIONS sequential big_endian'

  WRITE (nchout0,'(a,i10)')                                             &
      'FILEHEADER ',192        !!! The number depends on the file
                               !!! structure of sfcoutfl. See iolib3d.f

  WRITE (nchout0,'(a/a)')                                               &
      'UNDEF   -9.e+33','*'

  IF ( mapproj == 2 .OR. mapproj == -2) THEN
    WRITE (nchout0,'(a)')                                               &
        '* For lat-lon-lev display, umcomment the following 4 lines.'

    WRITE (nchout0,'(a,1x,i8,1x,i3,a,2f12.6,2i3,3f12.6,2f12.2)')        &
        'PDEF',nx,ny,' LCC',lat11,lon11,1,1,                            &
              trulat1,trulat2,trulon,xinc,yinc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'XDEF',nx, '  LINEAR  ',lonmin,loninc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'YDEF',ny, '  LINEAR  ',latmin,latinc
  ELSE
    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'XDEF',nx, '  LINEAR  ',xbgn,xinc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'YDEF',ny, '  LINEAR  ',ybgn,yinc
  END IF

  WRITE (nchout0,'(a,1x,i8,a,i1)')                                      &
      'ZDEF',1,'  LEVELS  ',0

  WRITE (nchout0,'(a,1x,i8,a,a,3x,i2.2,a/a)')                           &
      'TDEF',ntm,'  LINEAR  ',chrstr,tinc,dtunit,'*'

  WRITE (nchout0,'(a,1x,i3)')                                           &
      'VARS',varnum

  DO i = 1, varnum
    WRITE (nchout0,'(a8,1x,i3,2x,a,2x,a)')                              &
        varnam(i),varlev(i),varparam(i),vartit(i)
  END DO

  WRITE (nchout0,'(a)')                                                 &
      'ENDVARS'

  CLOSE (nchout0)
  CALL retunit(nchout0)

  RETURN
END SUBROUTINE sfccntl
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE SOILCNTL                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE soilcntl(nx,ny,nzsoil,zpsoil,soiloutfl,                     &
           zpsoilout,tsoilout,qsoilout,wcanpout,snowdout,x,y)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  ARPS soil data description file for GrADS display
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  10/05/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   Number of grid points in the soil
!
!    x        The x-coord. of the physical and computational grid.
!             Defined at u-point.
!    y        The y-coord. of the physical and computational grid.
!             Defined at v-point.
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny             ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  CHARACTER (LEN=*) :: soiloutfl ! Surface data file name

  INTEGER :: zpsoilout,tsoilout,qsoilout,wcanpout,snowdout
  INTEGER :: stypout

  REAL :: x(nx)                ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y(ny)                ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.

  REAL :: zpsoil(nx,ny,nzsoil) ! Depth of soil layers.

  REAL :: temxy1(nx,ny)        ! Temporary array
  REAL :: temxy2(nx,ny)        ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: varnumax
  PARAMETER ( varnumax = 100 )

  CHARACTER (LEN=256) :: soilctlfl
  CHARACTER (LEN=15) :: timestr,chrstr, chr1
  CHARACTER (LEN=8) :: varnam(varnumax)
  CHARACTER (LEN=60) :: vartit(varnumax)
  CHARACTER (LEN=10) :: varparam(varnumax)

  INTEGER :: varlev(varnumax)

  INTEGER :: varnum
  INTEGER :: nchout0

  CHARACTER (LEN=3) :: monnam(12)            ! Name of months
  DATA monnam/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                 &
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/

  CHARACTER (LEN=2) :: dtunit
  INTEGER :: ntm,tinc

  REAL :: lonmin,latmin,lonmax,latmax
  REAL :: xbgn,ybgn
  REAL :: xinc,yinc
  REAL :: lat11,lat12,lon11,lon12,lat21,lat22,lon21,lon22
  REAL :: latinc,loninc

  INTEGER :: ctllen, fnlen
  INTEGER :: i,j,k, is
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
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
!  Open the GrADS data control file
!
!-----------------------------------------------------------------------
!
  fnlen = LEN( soiloutfl )
  CALL strlnth( soiloutfl, fnlen )

  ctllen = fnlen + 4
  soilctlfl(1:ctllen) = soiloutfl(1:fnlen)//'.ctl'

  CALL fnversn( soilctlfl, ctllen )
  CALL getunit (nchout0)
  WRITE (6,'(a)') 'The GrADS data control file is '                     &
                    //soilctlfl(1:ctllen)

  OPEN (nchout0, FILE = soilctlfl(1:ctllen), STATUS = 'unknown')

  xbgn = 0.5 * (x(1) + x(2))
  ybgn = 0.5 * (y(1) + y(2))

  xinc = (x(2) - x(1))
  yinc = (y(2) - y(1))

  CALL xytoll(nx,ny,x,y,temxy1,temxy2)

  CALL xytoll(1,1,xbgn,ybgn,lat11,lon11)

  CALL a3dmax0(temxy1,1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,                    &
               latmax,latmin)
  CALL a3dmax0(temxy2,1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,                    &
               lonmax,lonmin)

  latinc = (latmax-latmin)/(ny-1)
  loninc = (lonmax-lonmin)/(nx-1)

  WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                                 &
           'latmin:latmax:latinc = ',                                   &
            latmin,':',latmax,':',latinc
  WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                                 &
           'lonmin:lonmax:loninc = ',                                   &
           lonmin,':',lonmax,':',loninc

  WRITE (timestr,'(i2.2,a,i2.2,a,i2.2,a3,i4.4)')                         &
      hour,':',minute,'Z',day,monnam(month),year

  ntm = 1
  tinc = 1
  dtunit = 'MN'

  DO i=1,varnumax
    varlev(i) = 0
    varparam(i) = '99'
  END DO

  varnum = 0
!
! 06/19/2002 Zuwen He
!
! Code changed due to new ou soil model.
!
! IF ( tsoilout+qsoilout+wcanpout /= 0 ) THEN
!   stypout = 1
! ELSE
!   stypout = 0
! END IF

  IF ( zpsoilout == 1 ) THEN
     varnum = varnum + 1
     varnam(varnum) = 'zpsoil'
     vartit(varnum) = 'Soil layer depth (m) '
     varlev(varnum) = nzsoil
  END IF


  IF ( tsoilout == 1 ) THEN
    DO is=0,nstyp
      WRITE (chrstr,'(i2.2)') is
      varnum = varnum + 1
      varnam(varnum) = 'tsoil'//chrstr(1:2)
      vartit(varnum) = 'Soil temperature (K) for soil type '//chrstr(1:2)
      if (is == 0) vartit(varnum) = 'Soil temperature (K)'
      varlev(varnum) = nzsoil
    END DO
  END IF

  IF ( qsoilout == 1 ) THEN
    DO is=0,nstyp
      WRITE (chrstr,'(i2.2)') is
      varnum = varnum + 1
      varnam(varnum) = 'qsoil'//chrstr(1:2)
      vartit(varnum) = 'Soil moisture (m**3/m**3) for soil type '//chrstr(1:2)
      if (is == 0) vartit(varnum) = 'Soil moisture (m**3/m**3)'
      varlev(varnum) = nzsoil
    END DO
  END IF

  IF ( wcanpout == 1 ) THEN
    DO is=0,nstyp
      WRITE (chrstr,'(i2.2)') is
      varnum = varnum + 1
      varnam(varnum) = 'wr'//chrstr(1:2)
      vartit(varnum) = 'Canopy water amount (m) for soil type '//chrstr(1:2)
        if (is == 0) vartit(varnum) = 'Canopy water amount (m)'
      varlev(varnum) = 0
    END DO
  END IF

  IF ( snowdout == 1 ) THEN
    varnum = varnum + 1
    varnam(varnum) = 'snowd'
    vartit(varnum) = 'Snow depth (m)'
    varlev(varnum) = 0
  END IF

! IF ( stypout == 1 ) THEN
!   DO is=1,nstyp
!     WRITE (chrstr,'(i2.2)') is
!     varnum = varnum + 1
!     varnam(varnum) = 'styp'//chrstr(1:2)
!     vartit(varnum) = 'Soil type '//chrstr(1:2)
!     varlev(varnum) = 0
!   END DO
! END IF

  WRITE (nchout0,'(a/a)')                                               &
      'TITLE   ARPS soil variable data control for '                    &
      //runname(1:lfnkey),'*'

  WRITE (nchout0,'(a,a)')                                               &
      'DSET    ',soiloutfl(1:fnlen)

!  PRINT *,'TEST #1 = fnlen = ',fnlen

  WRITE (nchout0,'(a)')                                                 &
      'OPTIONS sequential big_endian'

  WRITE (nchout0,'(a,i10)')                                             &
      'FILEHEADER ',236          !!! The number depends on the file
                                   !!! structure of soiloutfl. See iolib3d.f

  WRITE (nchout0,'(a/a)')                                               &
      'UNDEF   -9.e+33','*'

  IF ( mapproj == 2 .OR. mapproj == -2) THEN
    WRITE (nchout0,'(a)')                                               &
        '* For lat-lon-lev display, umcomment the following 4 lines.'

    WRITE (nchout0,'(a,1x,i8,1x,i3,a,2f12.6,2i3,3f12.6,2f12.2)')        &
        'PDEF',nx,ny,' LCC',lat11,lon11,1,1,                            &
              trulat1,trulat2,trulon,xinc,yinc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'XDEF',nx, '  LINEAR  ',lonmin,loninc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'YDEF',ny, '  LINEAR  ',latmin,latinc
  ELSE
    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'XDEF',nx, '  LINEAR  ',xbgn,xinc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'YDEF',ny, '  LINEAR  ',ybgn,yinc
  END IF

  WRITE (nchout0,'(a,1x,i8,a)')                                      &
      'ZDEF',nzsoil,'  LEVELS  '
!
! 06/19/2002 Zuwen He
!
! The soil model is upside down.
!
  WRITE (nchout0,'(8f10.2)')                                        &
    (zpsoil(1,1,k)-zpsoil(1,1,1),k=1,nzsoil)


  WRITE (nchout0,'(a,1x,i8,a,a,3x,i2.2,a/a)')                           &
      'TDEF',ntm,'  LINEAR  ',timestr,tinc,dtunit,'*'

  WRITE (nchout0,'(a,1x,i3)')                                           &
      'VARS',varnum

  DO i = 1, varnum
    WRITE (nchout0,'(a8,1x,i3,2x,a,2x,a)')                              &
        varnam(i),varlev(i),varparam(i),vartit(i)
  END DO

  WRITE (nchout0,'(a)')                                                 &
      'ENDVARS'

  CLOSE (nchout0)
  CALL retunit(nchout0)

  RETURN
END SUBROUTINE soilcntl
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE TRNCNTL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE trncntl(nx,ny, ternfn, x,y)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  ARPS terrain data description file for GrADS display
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  10/05/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    x        The x-coord. of the physical and computational grid.
!             Defined at u-point.
!    y        The y-coord. of the physical and computational grid.
!             Defined at v-point.
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny             ! Number of grid points in 3 directions

  CHARACTER (LEN=*) :: ternfn  ! Terrain data file name

  REAL :: x(nx)                ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y(ny)                ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.

  REAL :: temxy1(nx,ny)        ! Temporary array
  REAL :: temxy2(nx,ny)        ! Temporary array

  REAL :: xs(nx)        ! Temporary array
  REAL :: ys(ny)        ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: trnctlfl
  CHARACTER (LEN=15) :: chrstr, chr1

  INTEGER :: nchout0

  CHARACTER (LEN=3) :: monnam(12)            ! Name of months
  DATA monnam/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                 &
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/

  CHARACTER (LEN=2) :: dtunit
  INTEGER :: ntm,tinc

  REAL :: lonmin,latmin,lonmax,latmax
  REAL :: xbgn,ybgn
  REAL :: xinc,yinc
  REAL :: lat11,lat12,lon11,lon12,lat21,lat22,lon21,lon22
  REAL :: latinc,loninc

  INTEGER :: ctllen, fnlen
  INTEGER :: i,j, k, kk, is
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
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
!  Open the GrADS data control file
!
!-----------------------------------------------------------------------
!
  fnlen = LEN( ternfn )
  CALL strlnth( ternfn, fnlen )

  ctllen = fnlen + 4
  trnctlfl(1:ctllen) = ternfn(1:fnlen)//'.ctl'

  CALL fnversn( trnctlfl, ctllen )
  CALL getunit (nchout0)
  WRITE (6,'(a)') 'The GrADS data control file is '                     &
                    //trnctlfl(1:ctllen)

  OPEN (nchout0, FILE = trnctlfl(1:ctllen), STATUS = 'unknown')

  xbgn = 0.5 * (x(1) + x(2))
  ybgn = 0.5 * (y(1) + y(2))

  xinc = (x(2) - x(1))
  yinc = (y(2) - y(1))

  DO i = 1, nx-1
   xs(i) = (x(i)+x(i+1))*0.5
  END DO
  xs(nx) = x(nx)+xinc*0.5

  DO j = 1, ny-1
    ys(j) = (y(j)+y(j+1))*0.5
  END DO
  ys(ny) = y(ny)+yinc*0.5

  CALL xytoll(nx,ny,xs,ys,temxy1,temxy2)

  CALL xytoll(1,1,xbgn,ybgn,lat11,lon11)

  CALL a3dmax0(temxy1,1,nx,1,nx,1,ny,1,ny,1,1,1,1,                    &
               latmax,latmin)
  CALL a3dmax0(temxy2,1,nx,1,nx,1,ny,1,ny,1,1,1,1,                    &
               lonmax,lonmin)

  latinc = (latmax-latmin)/(ny-1)
  loninc = (lonmax-lonmin)/(nx-1)

  WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                                 &
           'latmin:latmax:latinc = ',                                   &
            latmin,':',latmax,':',latinc
  WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                                 &
           'lonmin:lonmax:loninc = ',                                   &
           lonmin,':',lonmax,':',loninc

  IF(month <= 0) month = 1
  WRITE (chrstr,'(i2.2,a,i2.2,a,i2.2,a3,i4.4)')                         &
      hour,':',minute,'Z',day,monnam(month),year

  ntm = 1
  tinc = 1
  dtunit = 'MN'

  WRITE (nchout0,'(a/a)')                                               &
      'TITLE   ARPS terrain data control for '                          &
      //runname(1:lfnkey),'*'

  WRITE (nchout0,'(a,a)')                                               &
      'DSET    ',ternfn(1:fnlen)

  WRITE (nchout0,'(a)')                                                 &
      'OPTIONS sequential big_endian'

  WRITE (nchout0,'(a,i10)')                                             &
      'FILEHEADER ',192        !!! The number depends on the file
                               !!! structure of ternfn. See iolib3d.f90

  WRITE (nchout0,'(a/a)')                                               &
      'UNDEF   -9.e+33','*'

  IF ( mapproj == 2 .OR. mapproj == -2) THEN
!    WRITE (nchout0,'(a)')                                               &
!        '* For lat-lon-lev display, umcomment the following 4 lines.'

    WRITE (nchout0,'(a,1x,i8,1x,i3,a,2f12.6,2i3,3f12.6,2f12.2)')        &
        'PDEF',nx,ny,' LCC',lat11,lon11,1,1,                            &
              trulat1,trulat2,trulon,xinc,yinc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'XDEF',nx, '  LINEAR  ',lonmin,loninc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'YDEF',ny, '  LINEAR  ',latmin,latinc

  ELSE IF (mapproj == 3) THEN  ! Mercator map projection
  
    CALL xytoll(1,1,xbgn+xinc,ybgn,lat21,lon21)

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'XDEF',nx, '  LINEAR  ',lon11,lon21-lon11

    WRITE (nchout0,'(a,1x,i8,a)') 'YDEF',ny, '  LEVELS'

    kk = 0
    DO j = 1,ny
          kk = kk+1
          WRITE(nchout0,'(F8.3,a)',ADVANCE='NO') temxy1(2,j)
          IF (kk == ny) THEN
            WRITE(nchout0,'(a)') ' '
          ELSE
            WRITE(nchout0,'(a)',ADVANCE='NO') ','
            IF (MOD(kk,10) == 0 ) WRITE(nchout0,'(a)') ' '
          END IF
    END DO

  ELSE
    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'XDEF',nx, '  LINEAR  ',xbgn,xinc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'YDEF',ny, '  LINEAR  ',ybgn,yinc
  END IF

  WRITE (nchout0,'(a,1x,i8,a,i1)')                                      &
      'ZDEF',1,'  LEVELS  ',0

  WRITE (nchout0,'(a,1x,i8,a,a,3x,i2.2,a/a)')                           &
      'TDEF',ntm,'  LINEAR  ',chrstr,tinc,dtunit,'*'

  WRITE (nchout0,'(a,1x,i3)')                                           &
      'VARS',1

  WRITE (nchout0,'(a)')                                                 &
      'trn      0   99   ARPS terrain (m)'

  WRITE (nchout0,'(a)')                                                 &
      'ENDVARS'

  CLOSE (nchout0)
  CALL retunit(nchout0)

  RETURN
END SUBROUTINE trncntl

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CHECKGRID3D               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE checkgrid3d(nx,ny,nz,nxin,nyin,nzin,                         &
           dx,dy,dz,dzmin,ctrlat,ctrlon,                                &
           strhopt,zrefsfc,dlayer1,dlayer2,zflat,strhtune,              &
           mapproj,trulat1,trulat2,trulon,sclfct,                       &
           dxin,dyin,dzin,dzminin,ctrlatin,ctrlonin,                    &
           strhoptin,zrefsfcin,dlayer1in,dlayer2in,zflatin,strhtunein,  &
           mapprojin,trulat1in,trulat2in,trulonin,sclfctin,ireturn)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check grid information to see if input data is consistent with
!  ARPS grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/24
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx
!  ny
!  nz
!  nxin
!  nyin
!  nzin
!  dx
!  dy
!  dz
!  dzmin
!  ctrlat
!  ctrlon
!  strhopt
!  zrefsfc
!  dlayer1
!  dlayer2
!  zflat
!  strhtune
!  mapproj
!  trulat1
!  trulat2
!  trulon
!  sclfct
!  dxin
!  dyin
!  dzin
!  dzminin
!  ctrlatin
!  ctrlonin
!  strhoptin
!  zrefsfcin
!  dlayer1in
!  dlayer2in
!  zflatin
!  strhtunein
!  mapprojin
!  trulat1in
!  trulat2in
!  trulonin
!  sclfctin
!
!  OUTPUT:
!
!  ireturn  Flag indicating if the grids are the same (0-okay,
!           1-significant differences)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nxin
  INTEGER :: nyin
  INTEGER :: nzin
  INTEGER :: nx,ny,nz
  REAL :: dx
  REAL :: dy
  REAL :: dz
  REAL :: dzmin
  REAL :: ctrlat
  REAL :: ctrlon
  INTEGER :: strhopt
  REAL :: zrefsfc
  REAL :: dlayer1
  REAL :: dlayer2
  REAL :: zflat
!  INTEGER :: strhtune
  REAL :: strhtune
  INTEGER :: mapproj
  REAL :: trulat1
  REAL :: trulat2
  REAL :: trulon
  REAL :: sclfct
  REAL :: dxin
  REAL :: dyin
  REAL :: dzin
  REAL :: dzminin
  REAL :: ctrlatin
  REAL :: ctrlonin
  INTEGER :: strhoptin
  REAL :: zrefsfcin
  REAL :: dlayer1in
  REAL :: dlayer2in
  REAL :: zflatin
!  INTEGER :: strhtunein
  REAL :: strhtunein
  INTEGER :: mapprojin
  REAL :: trulat1in
  REAL :: trulat2in
  REAL :: trulonin
  REAL :: sclfctin

  INTEGER :: ireturn

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  REAL :: eps
  PARAMETER (eps= 0.01)

!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Compare all the variables.  (Note that for 3-D data sets one ought
!  to also check hterain, which is not done here.)
!
!-----------------------------------------------------------------------

  ireturn = 0
  IF ((nxin /= nx) .OR. (nyin /= ny) .OR. (nzin /= nz)) THEN
    WRITE (6,'(/2(/5x,a),2(/5x,3(a,i9))/)')                             &
        'ERROR: nx, ny or nz in file does not match',                   &
        'the values in the program.',                                   &
        'In file,       nx=', nxin, ', ny=', nyin, ', nz=',nzin,        &
        'In executable, nx=', nx,   ', ny=', ny,   ', nz=',nz
    ireturn = 1
  END IF

  IF(ABS(dxin-dx) > eps .OR.  ABS(dyin-dy) > eps .OR.                   &
         ABS(dzin-dz) > eps .OR.                                        &
         ABS(dzminin-dzmin) > eps ) THEN
    WRITE(6,'(/2(/5x,a),4(/5x,2(a,f13.3))/)')                           &
        'ERROR: dx, dy, dz or dzmin in the file ',                      &
        'do not match those specified in the input file.',              &
        'In file,     dx=', dxin, ', dy=', dyin,                        &
        '             dz=', dzin, ', dzmin=', dzminin,                  &
        'In input,    dx=', dx,   ', dy=', dy,                          &
        '             dz=', dz, ', dzmin=', dzmin
    ireturn = 1
  END IF

  IF(ABS(ctrlatin-ctrlat) > eps .OR.  ABS(ctrlonin-ctrlon) > eps ) THEN
    WRITE(6,'(/2(/5x,a),2(/5x,2(a,f13.3))/)')                           &
        'ERROR: Central latitude and/or longitude of the grid ',        &
        'in the file do not match those specified in input file.',      &
        'In file,     ctrlat=',ctrlatin,', ctrlon=',ctrlonin,           &
        'In input,    ctrlat=',ctrlat  ,', ctrlon=',ctrlon
    ireturn = 1
  END IF

  IF(strhoptin /= strhopt .OR. ABS(zrefsfcin-zrefsfc) > eps .OR.        &
        ABS(dlayer1in-dlayer1) > eps .OR.                               &
        ABS(dlayer2in-dlayer2) > eps .OR.                               &
        ABS(zflatin-zflat) > eps .OR.                                   &
        ABS(strhtunein-strhtune) > eps .OR.                             &
        mapprojin /= mapproj .OR.                                       &
        ABS(trulat1in-trulat1) > eps .OR.                               &
        ABS(trulat2in-trulat2) > eps .OR.                               &
        ABS(trulonin-trulon ) > eps .OR.                                &
        ABS(sclfctin-sclfct ) > eps) THEN
    WRITE(6,'(/2(/5x,a),2(/5x,3(a,f13.3),2(a,i3),5(a,f13.3))/)')        &
         'ERROR: Map projection or other grid parameters do not ',      &
         'match those specified in input file.',                        &
         'In file,     trulat1=',trulat1in,', trulat2=',trulat2in,      &
         ', trulon=',trulonin,', mapproj=',mapprojin,                   &
         ', strhopt=',strhoptin,', zrefsfc=',zrefsfcin,                 &
         ', dlayer1=',dlayer1in,', dlayer2=',dlayer2in,                 &
         ', zflat=',zflatin,', strhtune=',strhtunein,                   &
         'In input,    trulat1=',trulat1  ,', trulat2=',trulat2,        &
         ', trulon=',trulon,   ', mapproj=',mapproj,                    &
         ', strhopt=',strhopt,', zrefsfc=',zrefsfc,                     &
         ', dlayer1=',dlayer1,', dlayer2=',dlayer2,                     &
         ', zflat=',zflat,', strhtune=',strhtune
    ireturn = 1
  END IF

  RETURN
END SUBROUTINE checkgrid3d

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CHECKGRID2D               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE checkgrid2d(nx,ny,nxin,nyin,                                 &
           dx,dy,ctrlat,ctrlon,                                         &
           mapproj,trulat1,trulat2,trulon,sclfct,                       &
           dxin,dyin,ctrlatin,ctrlonin,                                 &
           mapprojin,trulat1in,trulat2in,trulonin,sclfctin,ireturn)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check grid information to see if input data is consistent with
!  ARPS grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/03/24
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx
!  ny
!  nxin
!  nyin
!  dx
!  dy
!  ctrlat
!  ctrlon
!  mapproj
!  trulat1
!  trulat2
!  trulon
!  sclfct
!  dxin
!  dyin
!  ctrlatin
!  ctrlonin
!  mapprojin
!  trulat1in
!  trulat2in
!  trulonin
!  sclfctin
!
!  OUTPUT:
!
!  ireturn  Flag indicating if the grids are the same (0-okay,
!           1-significant differences)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nxin
  INTEGER :: nyin
  INTEGER :: nx,ny
  REAL :: dx
  REAL :: dy
  REAL :: ctrlat
  REAL :: ctrlon
  INTEGER :: mapproj
  REAL :: trulat1
  REAL :: trulat2
  REAL :: trulon
  REAL :: sclfct
  REAL :: dxin
  REAL :: dyin
  REAL :: ctrlatin
  REAL :: ctrlonin
  INTEGER :: mapprojin
  REAL :: trulat1in
  REAL :: trulat2in
  REAL :: trulonin
  REAL :: sclfctin

  INTEGER :: ireturn

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  REAL :: eps
  PARAMETER (eps= 0.1)

!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Compare all the variables.
!
!-----------------------------------------------------------------------

  ireturn = 0
  IF ((nxin /= nx) .OR. (nyin /= ny)) THEN
    WRITE (6,'(/2(/5x,a),2(/5x,(a,i9))/)')                              &
        'ERROR: nx and/or ny in file does not match',                   &
        'the values in the program.',                                   &
        'In file,       nx=', nxin, ', ny=', nyin,                      &
        'In executable, nx=', nx,   ', ny=', ny
    ireturn = 1
  END IF

  IF(ABS(dxin-dx) > eps .OR. ABS(dyin-dy) > eps) THEN
    WRITE(6,'(/2(/5x,a),2(/5x,2(a,f13.3))/)')                           &
        'ERROR: dx or dy in the file ',                                 &
        'do not match those specified in the input file.',              &
        'In file,     dx=', dxin, ', dy=', dyin,                        &
        'In input,    dx=', dx,   ', dy=', dy
    ireturn = 1
  END IF

  IF(ABS(ctrlatin-ctrlat) > eps .OR.  ABS(ctrlonin-ctrlon) > eps ) THEN
    WRITE(6,'(/2(/5x,a),2(/5x,2(a,f13.3))/)')                           &
        'ERROR: Central latitude and/or longitude of the grid ',        &
        'in the file do not match those specified in input file.',      &
        'In file,     ctrlat=',ctrlatin,', ctrlon=',ctrlonin,           &
        'In input,    ctrlat=',ctrlat  ,', ctrlon=',ctrlon
    ireturn = 1
  END IF

  IF(mapprojin /= mapproj .OR. ABS(trulat1in-trulat1) > eps .OR.        &
        ABS(trulat2in-trulat2) > eps .OR.                               &
        ABS(trulonin-trulon ) > eps .OR.                                &
        ABS(sclfctin-sclfct ) > eps) THEN
    WRITE(6,'(/2(/5x,a),2(/5x,3(a,f13.3),(a,i3,a,f13.3))/)')  &
         'ERROR: Map projection or other grid parameters do not ',      &
         'match those specified in input file.',                        &
         'In file,     trulat1=',trulat1in,', trulat2=',trulat2in,      &
         ', trulon=',trulonin,', mapproj=',mapprojin,' sclfct=',sclfctin, &
         'In input,    trulat1=',trulat1  ,', trulat2=',trulat2,        &
         ', trulon=',trulon,   ', mapproj=',mapproj,' sclfct=',sclfct
    ireturn = 1
  END IF

  RETURN
END SUBROUTINE checkgrid2d


!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE INITZTIME                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE initztime(ayear,amm,aday)


!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To calculate Zulu time based upon Mesonet observations read directly
!      for any given day.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: J.A. Brotzge
!  2001/08/23
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!
!    curtim    Computational time (seconds)
!
!  OUTPUT:
!
!    ztime     Zulu time, calculated (seconds)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  IMPLICIT NONE             ! Force explicit declarations

!
!  INPUT
!

!  INTEGER :: hour, minute, second   ! Computational time
!  REAL :: curtim                    ! Computational time (seconds)

!
!  OUTPUT
!

!  REAL :: ztime               !Zulu time (seconds)
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
  INTEGER :: dayt, initday
  CHARACTER :: ayear*4, amm*2, aday*2
  REAL :: modelinitime, day1time, day2time, day3time, day4time
  REAL :: day5time, day6time, day7time, day8time
  REAL :: day9time, day10time



!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  modelinitime = real(hour*3600 + minute*60 + second)
    day1time = 86400.0 - modelinitime
    day2time = day1time + 86400.0
    day3time = day2time + 86400.0
    day4time = day3time + 86400.0
    day5time = day4time + 86400.0
    day6time = day5time + 86400.0
    day7time = day6time + 86400.0
    day8time = day7time + 86400.0
    day9time = day8time + 86400.0
    day10time= day9time + 86400.0


    IF (curtim.lt.day1time) THEN
       dayt = day
      ztime = modelinitime + curtim
    ENDIF
    IF (curtim.ge.day1time.and.curtim.lt.day2time) THEN
       dayt = day + 1
      ztime = curtim - day1time
    ENDIF
    IF (curtim.ge.day2time.and.curtim.lt.day3time) THEN
       dayt = day + 2
      ztime = curtim - day2time
    ENDIF
    IF (curtim.ge.day3time.and.curtim.lt.day4time) THEN
       dayt = day + 3
      ztime = curtim - day3time
    ENDIF
    IF (curtim.ge.day4time.and.curtim.lt.day5time) THEN
       dayt = day + 4
      ztime = curtim - day4time
    ENDIF
    IF (curtim.ge.day5time.and.curtim.lt.day6time) THEN
       dayt = day + 5
      ztime = curtim - day5time
    ENDIF
    IF (curtim.ge.day6time.and.curtim.lt.day7time) THEN
       dayt = day + 6
      ztime = curtim - day6time
    ENDIF
    IF (curtim.ge.day7time.and.curtim.lt.day8time) THEN
       dayt = day + 7
      ztime = curtim - day7time
    ENDIF
    IF (curtim.ge.day8time.and.curtim.lt.day9time) THEN
       dayt = day + 8
      ztime = curtim - day8time
    ENDIF
    IF (curtim.ge.day9time.and.curtim.lt.day10time) THEN
       dayt = day + 9
      ztime = curtim - day9time
    ENDIF
    IF (curtim.ge.day10time) THEN
       dayt = day + 10
      ztime = curtim - day10time
    ENDIF

  IF (month.eq.2.and.day.gt.29) THEN
      month=month+1
      dayt=1
    ENDIF
      IF (month.eq.2.and.day.eq.29) THEN
        IF (year.ne.2000) THEN
          month=month+1
          dayt=1
        ENDIF
      ENDIF
    IF (dayt.eq.31) THEN
      IF (month.eq.4.or.month.eq.6) THEN
        month=month+1
        dayt=1
      ENDIF
      IF (month.eq.9.or.month.eq.11) THEN
        month=month+1
        dayt=1
      ENDIF
    ENDIF
    IF (dayt.eq.32) THEN
      month=month+1
      dayt=1
    ENDIF

   write(ayear,'(i4)') year
   write(amm,'(i2)') month
   write(aday,'(i2)') dayt
   IF (month.lt.10) amm(1:1)='0'
   IF (dayt.lt.10) aday(1:1)='0'

  RETURN
END SUBROUTINE initztime


!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READJSOIL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readjsoil(nx,ny,nzsoil,numbsty,ayear,amm,aday,arpstime,  &
           zpsoil,tsoil,qsoil,wetcanp,snowdpth)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in the mesonet radt'n data qc'd and processed by Jerry Brotzge.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber and Jerry Brotzge
!  8/28/2001.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   Number of soil layers
!    numbsty  Number of soil types within grid box
!
!  OUTPUT:
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m^3/m^3)
!    wetcanp  Canopy soil moisture (unitless)
!    snowdpth Snow depth (m)
!    soiltyp  Soil type (unitless)
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  IMPLICIT NONE             ! Force explicit declarations

!
!  INPUT
!

  INTEGER :: nx, ny            ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil layers

  INTEGER :: numbsty             ! Number of soil types within grid box

  INTEGER, PARAMETER :: nmeso=48, mmeso=288
  INTEGER, PARAMETER :: nzsoilext=5
  INTEGER, PARAMETER :: nzsoil_ext=5


    REAL :: t105 (nmeso)
    REAL :: t125 (nmeso)
    REAL :: t160 (nmeso)
    REAL :: t175 (nmeso)
    REAL :: smo05 (nmeso)
    REAL :: smo25 (nmeso)
    REAL :: smo60 (nmeso)
    REAL :: smo75 (nmeso)
    REAL :: wrxx (nmeso)
    REAL :: qvsf (nmeso)
    REAL :: snow (nmeso)

  REAL :: rnet  (mmeso)     ! Mesonet CNR1 rnet (W/m**2)
  REAL :: hflx  (mmeso)     ! Mesonet sonic sensible heat flux (W/m**2)
  REAL :: lflx  (mmeso)     ! Mesonet sonic latent heat flux (W/m**2)
  REAL :: gflx  (mmeso)     ! Mesonet ground heat flux (W/m**2)
  REAL :: l2fx  (mmeso)     ! Mesonet residual latent heat flux (W/m**2)
  REAL :: irtx  (mmeso)     ! Mesonet rainfall rate (kg/m^2/s)
  REAL :: irthour (mmeso)
  REAL :: prof  (mmeso)

  CHARACTER(4) :: temp1 (nmeso)  ! Mesonet station name
  INTEGER :: stnm1 (nmeso)      ! Mesonet station number
  INTEGER :: mesotime(nmeso)    ! Mesonet time (minutes then seconds)
  CHARACTER(4) :: ftemp1 (mmeso)  ! Mesonet station name
  INTEGER :: fstnm1 (mmeso)      ! Mesonet station number
  INTEGER :: fmesotime(mmeso)    ! Mesonet time (minutes then seconds)
  REAL :: mesot(nmeso)          ! Mesonet time (seconds)

!
!  OUTPUT
!

    REAL :: zpsoil(nx,ny,nzsoil)        ! Depth of soil (m)
    REAL :: tsoil(nx,ny,nzsoil,numbsty) ! Soil temperature (K)
    REAL :: qsoil(nx,ny,nzsoil,numbsty) ! Soil moisture (m3/m3)

    REAL :: wetcanp(nx,ny,numbsty) ! Canopy moisture
    REAL :: snowdpth(nx,ny)      ! Snow depth

    INTEGER :: soiltyp(nx,ny,numbsty)   !Soil type

!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,n,nt
  CHARACTER :: ayear*4, amm*2, aday*2

  CHARACTER :: atemp*4, btemp*2, ctemp*2, dtemp*3
  CHARACTER :: astid*4, astnm*4, aime*4

  CHARACTER :: a125*4, a160*4, a175*4, atm05*4, atm25*4
  CHARACTER :: atm60*4, atm75*4, awrxx*4, aqvsf*4, asnow*4
  CHARACTER :: a105*4
  CHARACTER :: arnet*4, ahflx*4, alflx*4, al2fx*4
  CHARACTER :: agflx*4, airtx*4
  CHARACTER :: bstid*4, bstnm*4, bime*4
  CHARACTER :: etemp*4, ftemp*4, gtemp*4, htemp*4

  INTEGER :: stnm, cmm, cdd
  INTEGER :: fmm, fdd, filelength

  REAL :: tema, temb, convertp
  REAL :: arpstime

  REAL :: zpsoil_ext(nx,ny,nzsoilext)
  REAL :: tsoil_ext(nx,ny,nzsoilext), qsoil_ext(nx,ny,nzsoilext)

!  DATA zpsoil_ext/0.0,-0.05,-0.25,-0.60,-0.75/

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'indtflg.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

    filelength = len_trim(sitesoil)

!    nzsoil_ext = 5

     OPEN(unit=60,file=sitesoil(1:filelength)//                      &
        ayear//amm//aday//'.mts',status='old',form='formatted')

        READ(60,311) dtemp
        READ(60,313) btemp, atemp, cmm, cdd, ctemp, ctemp, ctemp
        READ(60,301) astid, astnm, aime, a105, atm05,                 &
             & a125, atm25, a160, atm60, a175, atm75,                 &
             & awrxx, aqvsf, asnow


 301   FORMAT (a4,13(1x,a4))
 302   FORMAT (1x,a4,4x,i3,4x,i5,11(1x,f9.4))
 311   FORMAT (6x,a3)
 313   FORMAT (7x,a2,1x,a4,1x,i2,1x,i2,1x,a2,1x,a2,1x,a2)


    filelength = len_trim(siteflux)

     OPEN(unit=20,file=siteflux(1:filelength)//                       &
  &     ayear//amm//aday//'.mts',status='old',form='formatted')

        read(20,511) htemp
        read(20,513) ftemp, etemp, fmm, fdd, gtemp, gtemp, gtemp
        read(20,501) bstid, bstnm, bime, arnet, ahflx, alflx,       &
             agflx, al2fx, airtx

 501   FORMAT (a4,8(1x,a4))
 502   format (1x,a4,4x,i3,5x,i4,7(1x,f9.4))
 511   FORMAT (6x,a3)
 513   FORMAT (7x,a2,1x,a4,1x,i2,1x,i2,1x,a2,1x,a2,1x,a2)

     nt = 0

      DO n=1,288   ! perform the file read....

        read(20,502) ftemp1(n), fstnm1(n), fmesotime(n), &
     &        rnet(n), hflx(n), lflx(n),            &
     &        gflx(n), l2fx(n), irtx(n), prof(n)

     IF (n.eq.1.or.mod(n,6).eq.0.0) THEN
       nt = nt + 1
       irthour(nt) = irtx(n) + 273.15    !Convert C to K
     ENDIF

     ENDDO        !end of file read for Flux files

      DO j=1,ny
      DO i=1,nx
        zpsoil_ext(i,j,1) =  0.00
        zpsoil_ext(i,j,2) = -0.05
        zpsoil_ext(i,j,3) = -0.25
        zpsoil_ext(i,j,4) = -0.60
        zpsoil_ext(i,j,5) = -0.75
      END DO
      END DO


      DO n=1,48   ! perform the file read ....

        READ(60,302) temp1(n),stnm1(n),mesotime(n),              &
          t105(n),smo05(n),t125(n),smo25(n),                     &
          t160(n),smo60(n),t175(n),smo75(n),                     &
          wrxx(n),qvsf(n),snow(n)

        t105(n) = t105(n) + 273.15
        t125(n) = t125(n) + 273.15
        t160(n) = t160(n) + 273.15
        t175(n) = t175(n) + 273.15

        mesotime(n) = mesotime(n) * 60
        mesot(n) = real(mesotime(n))

  IF (arpstime.le.mesot(n)) THEN

        IF(arpstime.eq.mesot(n)) THEN   ! we have an exact match

        DO j=1,ny
        DO i=1,nx
          tsoil_ext(i,j,1) = irthour(n)
          tsoil_ext(i,j,2) = t105(n)
          tsoil_ext(i,j,3) = t125(n)
          tsoil_ext(i,j,4) = t160(n)
          tsoil_ext(i,j,5) = t175(n)

          qsoil_ext(i,j,1) = smo05(n)
          qsoil_ext(i,j,2) = smo05(n)
          qsoil_ext(i,j,3) = smo25(n)
          qsoil_ext(i,j,4) = smo60(n)
          qsoil_ext(i,j,5) = smo75(n)
        END DO
        END DO

        CALL initsoil(nx,ny,nzsoil,nzsoil_ext,numbsty,zpsoil,  &
                 zpsoil_ext,tsoil,tsoil_ext,qsoil,qsoil_ext)

          DO j = 1,ny
          DO i = 1,nx
            wetcanp(i,j,1) = wrxx(n)
            snowdpth(i,j) = snow(n)
          END DO
          END DO

        ELSE IF(arpstime.gt.mesot(n-1).and.arpstime.le.mesot(n)) THEN

          tema = (arpstime-mesot(n-1))/1800.0
          temb = 1.0-tema

          DO j = 1,ny
          DO i = 1,nx

            tsoil_ext(i,j,1) = temb*irthour(n-1)+tema*irthour(n)
            tsoil_ext(i,j,2) = temb*t105(n-1)+tema*t105(n)
            tsoil_ext(i,j,3) = temb*t125(n-1)+tema*t125(n)
            tsoil_ext(i,j,4) = temb*t160(n-1)+tema*t160(n)
            tsoil_ext(i,j,5) = temb*t175(n-1)+tema*t175(n)

            qsoil_ext(i,j,1) = temb*smo05(n-1)+tema*smo05(n)
            qsoil_ext(i,j,2) = temb*smo05(n-1)+tema*smo05(n)
            qsoil_ext(i,j,3) = temb*smo25(n-1)+tema*smo25(n)
            qsoil_ext(i,j,4) = temb*smo60(n-1)+tema*smo60(n)
            qsoil_ext(i,j,5) = temb*smo75(n-1)+tema*smo75(n)

            wetcanp(i,j,1) = temb*wrxx(n-1)+tema*wrxx(n)
            snowdpth(i,j) = temb*snow(n-1)+tema*snow(n)
          END DO
          END DO

        CALL initsoil(nx,ny,nzsoil,nzsoil_ext,numbsty,zpsoil,  &
                 zpsoil_ext,tsoil,tsoil_ext,qsoil,qsoil_ext)

        END IF   !  end of the time if loop
    END IF       !  end of the time if loop for matching arpstime and mesot

     ENDDO                                   !Time loop

        close (20)
        close (60)
        sfcin = 1
        return

  RETURN
 END SUBROUTINE readjsoil


!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READJVEG                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readjveg(nx,ny,numbsty,ayear,amm,aday,arpstime,         &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi)

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in the mesonet radt'n data qc'd and processed by Jerry Brotzge.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber and Jerry Brotzge
!  8/28/2001.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    numbsty  Number of soil types within grid box
!
!  OUTPUT:
!
!    soiltyp     Soil type
!    stypfrct    Surface type fraction
!    vegtyp      Vegetation type (unitless)
!    lai         Leaf area index ()
!    roufns      Roughness (unitless)
!    veg         Vegetation fraction (unitless)
!    ndvi        Normalized Diff Veg Index (unitless)
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!-----------------------------------------------------------------------
!
!

  IMPLICIT NONE             ! Force explicit declarations

!
!  INPUT
!

  INTEGER :: nx, ny            ! Number of grid points in 3 directions
  INTEGER :: numbsty           ! Number of soil types within grid box
  INTEGER :: nst               ! Number of soil types

  INTEGER, PARAMETER :: nmeso = 48

    REAL :: veg2 (nmeso)
    REAL :: veg4 (nmeso)
    REAL :: veg5 (nmeso)
    REAL :: veg6 (nmeso)
    REAL :: veg7 (nmeso)

    INTEGER :: veg1(nmeso)
    INTEGER :: veg3(nmeso)

  CHARACTER(4) :: temp1 (nmeso) ! Mesonet station name
  INTEGER :: stnm1 (nmeso)  ! Mesonet station number
  INTEGER :: mesotime(nmeso) !Mesonet time (minutes then seconds)
  REAL :: mesot(nmeso)   ! Mesonet time (seconds)

!
!  OUTPUT
!

    REAL :: stypfrct(nx,ny,numbsty)    !Surface type fraction
    REAL :: lai(nx,ny)               !Leaf Area Index
    REAL :: roufns(nx,ny)            !Surface roughness
    REAL :: veg(nx,ny)               !Vegetation
    REAL :: ndvi(nx,ny)              !Normalized Diff Veg Index

    INTEGER :: soiltyp(nx,ny,numbsty)   !Soil type
    INTEGER :: vegtyp(nx,ny)          !Veg type
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,n
  CHARACTER :: ayear*4, amm*2, aday*2

  CHARACTER :: atemp*4, btemp*2, ctemp*2, dtemp*3
  CHARACTER :: astid*4, astnm*4, aime*4

  CHARACTER :: astyp*4, asfct*4, avtyp*4, alaix*4, arfns*4
  CHARACTER :: avfct*4, andvi*4

  INTEGER :: stnm, cmm, cdd, filelength

  REAL :: tema, temb, convertp
  REAL :: arpstime
!
!-----------------------------------------------------------------------
!
!  Include files:
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

    filelength = len_trim(siteveg)

     OPEN(unit=70,file=siteveg(1:filelength)//                   &
  &    ayear//amm//aday//'.mts',status='old',form='formatted')

!    loop that reads Jerry's mesonet data file (qc'd)
!    note that the ARPS time must match the mesonet datafile time.

        READ(70,711) dtemp
        READ(70,713) btemp, atemp, cmm, cdd, ctemp, ctemp, ctemp
        READ(70,701) astid, astnm, aime, astyp, asfct,             &
                     avtyp, alaix, arfns, avfct, andvi


 701   FORMAT (a4,9(1x,a4))
 702   FORMAT (1x,a4,4x,i3,5x,i4,1x,i2,1x,f4.2,1x,i2,              &
           4(1x,f4.2))
 711   FORMAT (6x,a3)
 713   FORMAT (7x,a2,1x,a4,1x,a2,1x,a2,1x,a2,1x,a2,1x,a2)

      DO n=1,48   ! perform the file read....

        READ(70,702) temp1(n),stnm1(n),mesotime(n),              &
          veg1(n),veg2(n),veg3(n),veg4(n),                       &
          veg5(n),veg6(n),veg7(n)

        mesotime(n) = mesotime(n) * 60
        mesot(n) = real(mesotime(n))

  IF (arpstime.le.mesot(n)) THEN

        IF(arpstime.eq.mesot(n)) THEN   ! we have an exact match

          DO j = 1,ny
          DO i = 1,nx
            soiltyp(i,j,1) = veg1(n)
            stypfrct(i,j,1) = veg2(n)
            vegtyp(i,j) = veg3(n)
            lai(i,j) = veg4(n)
            roufns(i,j) = veg5(n)
            veg(i,j) = veg6(n)
            ndvi(i,j) = veg7(n)

          END DO
          END DO

        ELSE IF(arpstime.gt.mesot(n-1).and.arpstime.le.mesot(n))THEN

          tema = (arpstime-mesot(n-1))/1800.0
          temb = 1.0-tema

          DO j = 1,ny
          DO i = 1,nx
            soiltyp(i,j,1) = veg1(n)
            stypfrct(i,j,1) = veg2(n)
            vegtyp(i,j) = veg3(n)
            lai(i,j) = veg4(n)
            roufns(i,j) = veg5(n)
            veg(i,j) = veg6(n)
            ndvi(i,j) = veg7(n)

          END DO
          END DO

        END IF   !  end of the time if loop
    END IF       !  end of the time if loop for matching arpstime and mesot


      ENDDO                                  !Time loop

        close (70)
        return

  RETURN
END SUBROUTINE readjveg

!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE INITSOIL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initsoil (nx,ny,nzsoil,nzsoil_ext, nstyp,zpsoil,  &
         zpsoil_ext,tsoil,tsoil_ext,qsoil,qsoil_ext)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take soil data from external file and interpolate in
!  the vertical to form the ARPS 2D data set profile
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Jerry Brotzge
!  05/15/2002
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!             for the ARPS grid
!    ny       Number of grid points in the y-direction (north/south)
!             for the ARPS grid
!    nzsoil   Number of grid points in the soil
!
!    nzsoil_ext   Number of grid points in the soil external file
!
!    nstyp    Number of soil types
!
!    zpsoil   Physical depth of soil layers (m)
!
!    zpsoil_ext Physical depth of external model soil layers (m)
!
!    tsoil    Soil temperature (K)
!
!    qsoil    Soil moisture (m**3/m**3)
!
!  OUTPUT:
!
!    tsoil    Soil temperature profile (K)
!
!    qsoil    Soil moisture profile (m**3/m**3)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!
!-----------------------------------------------------------------------
!
!  Input/output variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kk

  INTEGER :: nx,ny                 ! ARPS grid dimensions
  INTEGER :: nzsoil                ! ARPS soil levels
  INTEGER :: nzsoil_ext            ! External file soil levels
  INTEGER :: nstyp                 ! Number of soil types
!
  REAL :: zpsoil(nx,ny,nzsoil)     ! Physical depth of ARPS soil layers
  REAL :: zpsoil_ext(nx,ny,nzsoil_ext) ! Physical depth of ext. file soil layers

  REAL :: tsoil(nx,ny,nzsoil)      ! Soil temperature (K)
  REAL :: qsoil(nx,ny,nzsoil)      ! Soil moisture (m**3/m**3)

  REAL :: tsoil_ext(nx,ny,nzsoil_ext) ! External file soil temperature (K)
  REAL :: qsoil_ext(nx,ny,nzsoil_ext) ! External file soil moisture (m**3/m**3)


!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: index
  REAL :: dampdepth
  REAL :: tsoil15
  REAL :: qsoil15
  REAL :: dampdz
  REAL :: depthdz(nx,ny,nzsoil)
  REAL :: depthdzext

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'grid.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!-----------------------------------------------------------------------
!
!    Vertical interpolation
!
!-----------------------------------------------------------------------
!
!  NOTE: Index INCREASES WITH DEPTH. k=1 at the top; nzsoil at the bottom BC.
!
!----------------------------------------------------------------
!     Set top boundary condition at level = 1
!-------------------------------------------------------------

  DO i=1,nx
    DO j=1,ny

      tsoil (i,j,1) = tsoil_ext(i,j,1)
      qsoil (i,j,1) = qsoil_ext(i,j,1)

      DO kk=1,nzsoil_ext
        zpsoil_ext(i,j,kk) = zpsoil_ext(i,j,kk) + zrefsfc
      END DO

    END DO
  END DO

!----------------------------------------------------------------------
!   Initialize soil temp and moisture profiles
!----------------------------------------------------------------------

! Note: All soil depths are negative (downward)

  DO k=2,nzsoil
    DO j=1,ny
      DO i=1,nx
        dampdepth = zpsoil(i,j,1) - 0.15  !Typical damping depth = 15 cm

        dampdz = zpsoil(i,j,k) - dampdepth
        depthdz(i,j,k) = zrefsfc - zpsoil(i,j,k)

        tsoil15 = 0.0
        qsoil15 = 0.0
        index = 0

        DO kk=2,nzsoil_ext

        depthdzext  = zrefsfc - zpsoil_ext(i,j,kk)

           IF (((zpsoil_ext(i,j,kk-1) >= dampdepth).AND.(zpsoil_ext(i,j,kk) <= &
                        dampdepth)).AND.(index < 1)) THEN

              tsoil15 = tsoil_ext(i,j,kk-1)+((tsoil_ext(i,j,kk)-   &
                        tsoil_ext(i,j,kk-1))* ((dampdepth - zpsoil_ext(i,j,kk))/ &
                        (zpsoil_ext(i,j,kk-1) - zpsoil_ext(i,j,kk))) )

              qsoil15 = qsoil_ext(i,j,kk-1)+((qsoil_ext(i,j,kk)-   &
                        qsoil_ext(i,j,kk-1))* ((dampdepth - zpsoil_ext(i,j,kk))/ &
                        (zpsoil_ext(i,j,kk-1) - zpsoil_ext(i,j,kk))) )

              index = index + 1
            ENDIF

         ENDDO

         DO kk=2,nzsoil_ext

!----------------------------------------------------------------------
!   Initialize levels equal to obs when levels are equal.
!----------------------------------------------------------------------

   IF (zpsoil(i,j,k) == zpsoil_ext(i,j,kk)) THEN

            tsoil(i,j,k) = tsoil_ext(i,j,kk)
            qsoil(i,j,k) = qsoil_ext(i,j,kk)

!----------------------------------------------------------------------
!   Exponential fit to initialization above damping depth (15 cm)
!----------------------------------------------------------------------

   ELSE IF (zpsoil(i,j,k) >= dampdepth) THEN  !Soil level within top 15 cm

          IF ((zpsoil(i,j,k) > zpsoil_ext(i,j,kk)).AND.(zpsoil_ext(i,j,kk) >=  &
              dampdepth)) THEN

                 tsoil(i,j,k) = tsoil_ext(i,j,kk)+((tsoil(i,j,k-1)-             &
                       tsoil_ext(i,j,kk))*EXP(- (depthdz(i,j,k)/depthdzext)) )

                 qsoil(i,j,k) = qsoil_ext(i,j,kk)+((qsoil(i,j,k-1)-             &
                       qsoil_ext(i,j,kk))*EXP(- (depthdz(i,j,k)/depthdzext)) )


          ELSE IF ((zpsoil(i,j,k) > zpsoil_ext(i,j,kk)).AND.(zpsoil_ext(i,j,kk) <  &
              dampdepth)) THEN   !(Est's linear fit upward to 15cm from tsoil_ext)


                tsoil(i,j,k)=tsoil15 + (tsoil(i,j,k-1)-tsoil15)* &
                        EXP( - (depthdz(i,j,k)/dampdz) )

                qsoil(i,j,k)=qsoil15 + (qsoil(i,j,k-1)-qsoil15)* &
                        EXP( - (depthdz(i,j,k)/dampdz) )


          ELSE IF ((zpsoil(i,j,k) < zpsoil_ext(i,j,kk-1)).AND.(zpsoil(i,j,k) >  &
              zpsoil_ext(i,j,kk))) THEN

                tsoil(i,j,k)=tsoil_ext(i,j,kk)+(tsoil_ext(i,j,kk-1) -             &
                   tsoil_ext(i,j,kk)) * EXP( - (depthdz(i,j,k)/depthdzext) )

                qsoil(i,j,k)=qsoil_ext(i,j,kk)+(qsoil_ext(i,j,kk-1) -             &
                   qsoil_ext(i,j,kk)) * EXP( - (depthdz(i,j,k)/depthdzext) )

            END IF


!----------------------------------------------------------------------
!   Linear fit between initialized levels below damping depth (15 cm)
!----------------------------------------------------------------------

        ELSE IF (zpsoil(i,j,k) < dampdepth) THEN  !Soil level below top 15 cm

          IF (zpsoil(i,j,k) < zpsoil_ext(i,j,kk-1) .AND.   &
                   zpsoil(i,j,k) > zpsoil_ext(i,j,kk)) THEN

            tsoil(i,j,k) = tsoil_ext(i,j,kk)+((tsoil_ext(i,j,kk-1)-   &
               tsoil_ext(i,j,kk)) * (zpsoil(i,j,k) - zpsoil_ext(i,j,kk))/ &
               (zpsoil_ext(i,j,kk-1) - zpsoil_ext(i,j,kk)) )

            qsoil(i,j,k) = qsoil_ext(i,j,kk)+((qsoil_ext(i,j,kk-1)-   &
               qsoil_ext(i,j,kk)) * (zpsoil(i,j,k) - zpsoil_ext(i,j,kk))/  &
               (zpsoil_ext(i,j,kk-1) - zpsoil_ext(i,j,kk)) )


          ELSE IF (zpsoil(i,j,k) > zpsoil_ext(i,j,2)) THEN

            tsoil(i,j,k) = tsoil_ext(i,j,kk)+((tsoil15 - tsoil_ext(i,j,kk)) * &
                ((zpsoil(i,j,k) - zpsoil_ext(i,j,kk))/(dampdepth - &
                  zpsoil_ext(i,j,kk))) )

            qsoil(i,j,k) = qsoil_ext(i,j,kk)+((qsoil15 - qsoil_ext(i,j,kk)) * &
                ((zpsoil(i,j,k) - zpsoil_ext(i,j,kk))/(dampdepth - &
                  zpsoil_ext(i,j,kk))) )

          ELSE IF (zpsoil(i,j,k) < zpsoil_ext(i,j,nzsoil_ext)) THEN

              tsoil(i,j,k) = tsoil_ext(i,j,kk)
              qsoil(i,j,k) = qsoil_ext(i,j,kk)

          END IF
      END IF        !Damping depth if/then

       END DO

       END DO
     END DO
   END DO

  RETURN
END SUBROUTINE initsoil

