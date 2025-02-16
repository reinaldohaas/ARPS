!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDNMCGRB                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdnmcgrb(nx_ext,ny_ext,nz_ext,                               &
                    gribfile,grbflen, gribtime,                         &
                    gridesc, iproj_grb, gthin,                          &
                    ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,      &
                    pi_grb,pj_grb,ipole, di_grb,dj_grb,                 &
                    latsw,lonsw, latne,lonne,                           &
                    latrot,lonrot,angrot,                               &
                    latstr,lonstr,facstr,                               &
                    lattru1,lattru2,lontrue,                            &
                    scanmode, iscan,jscan,kscan,                        &
                    ires,iearth,icomp,                                  &
                    jpenta,kpenta,mpenta,ispect,icoeff,                 &
                    xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                &
                    rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg, iret)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is to read and select given variables with GRIB
!  grid ID and level type from NMC GRIB files.
!
!  The decoder of GRIB is from NMC, ftp:nic.fb4.noaa.gov/pub
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/02/1996 Initialized.
!
!  MODIFICATIONS:
!
!  09/24/1998 (Dennis A. Moon, Chief Scientist, SSESCO)
!  Added sorting in vertical levels to assure that the level changes
!  monotonically decreasing.
!
!  1999/08/03 (Gene Bassett)
!  Excluded soil levels (grid types 111 & 112) from vertical level sorting.
!
!  2000/01/05 (Eric Kemp)
!  Fixed problem with century change over (ipds(21) runs from 1 to 100).
!
!  2002/09/04 (Kevin W. Thomas)
!  ETA soil temp and moisture data were never sorted, as they had flags
!  indicating "sigma vertical coordinates".  When CONDUIT (LDM) data is
!  the source, the data order isn't guaranteed, so it needs to be sorted.
!
!  Duplicated data reports that sometimes occur in CONDUIT data cause fatal
!  errors.  Duplicates are ignored.
!
!  2004/08/23 (Kevin W. Thomas)
!  Variable "grbunit" is now an allocatable character string.  It is
!  guaranteed to be the right size to hold a C file pointer.  This will
!  end the "segmentation violation" problems.
!
!  2006/10/19 (Kevin W. Thomas)
!  Avoid skipping record if PDS(16)=1.  This shows up in all 0 hour
!  records when NCEP moved the data from one system to another.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    gribfile      GRIB file name
!    grbflen       Length of GRIB file name
!    gribtime      Time of GRIB data
!
!  OUTPUT:
!
!    iproj_grb     Map projection number of GRIB data
!    trlon_grb     True longitude of GRIB data (degrees E)
!    latnot_grb(2) True latitude(s) of GRIB data (degrees N)
!    swlat_grb     Latitude  of first grid point at southwest conner
!    swlon_grb     Longitude of first grid point at southwest conner
!
!    dx_grb        x-direction grid length
!    dy_grb        y-direction grid length
!
!    iret          Return flag
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'gribcst.inc'      ! General GRIB definitions
  INCLUDE 'mp.inc'

  INTEGER :: nx_ext,ny_ext,nz_ext

  INTEGER :: grbflen, grbflen_new
  CHARACTER (LEN=*) :: gribfile
  CHARACTER (LEN=*) :: gribtime

  CHARACTER (LEN=256) :: gribfile_new
  CHARACTER (LEN=1)   :: opnmod

  INTEGER, INTENT(IN) :: lvldbg
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: gdsflg       ! GDS indicator

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth

  INTEGER :: iret         ! Return flag

!
!-----------------------------------------------------------------------
!
!  Temporary arrays to read GRIB file
!
!-----------------------------------------------------------------------
!
  INTEGER :: nrecs        ! number of records

  INTEGER :: fsize        ! Size of GRIB file
  INTEGER :: stdio_size   ! Size of "struct file" (C code in gribio_c.c).
  INTEGER :: istatus      ! Return status
  CHARACTER(len=1), allocatable :: grbunit(:) ! C structure, size stdio_size

  !wdt RLC
  LOGICAL (KIND=1) :: ibms(nx_ext*ny_ext)    ! BMS logical array for data bit map

  REAL    :: rcdata(nx_ext*ny_ext)           ! temporary data array
  REAL    :: var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt)        ! GRIB 2-D variables
  REAL    :: var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt) ! GRIB 3-D variables
  INTEGER :: var_lev3d(nz_ext,n3dvs,n3dlvt)
                                       ! Levels (hybrid) for each 3-D variable
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,l,m,n, ij, ki, kj, ilev
  INTEGER :: iyr,imo,iday,ihr,fhr, myr

  INTEGER :: ileft, iright, iincr
  INTEGER :: jleft, jright, jincr

  INTEGER :: chklev, mlvtp
  INTEGER :: iendn,itypec,wdlen, nwrd

  INTEGER :: itmp            ! another return flag

  LOGICAL :: grdflg           ! flag for grid information
  LOGICAL :: fexist

  INTEGER :: mismatch

  LOGICAL :: verbose = .FALSE.

  CHARACTER(LEN=11), PARAMETER :: cmd     = 'rdnmcgrib: '
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (lvldbg > 90) verbose = .TRUE.
  mismatch = 0
!
!-----------------------------------------------------------------------
!
!  Find the type of endian, type of character set and length of
!  machine word
!
!  iendn     -  Integer for big-endian or little-endian
!               = 0   big-endian
!               = 1   little-endian
!               = 2   cannot compute
!  itypec    -  Integer for type of character set
!               = 0   ascii  character set
!               = 1   ebcdic character set
!               = 2   not ascii or ebcdic
!  wdlen     -  Integer for words size of computer in bytes
!               = 4   for 32 bit (4 bytes) computers
!               = 8   for 64 bit (8 bytes) computers
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    CALL findvartype(iendn,itypec,wdlen)
    IF ( iendn == 2 ) THEN
      WRITE(6,'(a)') 'Unknown endian type. Stopped in rdnmcgrb'
      CALL arpsstop("endian error",1)
    END IF

    IF ( itypec == 2 ) THEN
      WRITE(6,'(a)') 'Unknown character set. Stopped in rdnmcgrb'
      CALL arpsstop("endian error",1)
    END IF
!
!-----------------------------------------------------------------------
!
!   initialize the var_lev3d array so we can sort the records by
!   level values as we read them. We can't assume monotonic ordering
!   of the levels
!
!-----------------------------------------------------------------------
!
    DO k=1, nz_ext
      DO n=1, n3dvs
        DO m=1, n3dlvt
          var_lev3d(k,n,m)= -999999
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Open the GRIB file
!
!-----------------------------------------------------------------------
!
    itmp = LEN_TRIM(gribtime)
    IF (itmp == 13) THEN
      READ (gribtime,'(i4,3i2,1x,i2)') iyr,imo,iday,ihr,fhr
    ELSE IF (itmp == 14) THEN
      READ (gribtime,'(i4,3i2,1x,i3)') iyr,imo,iday,ihr,fhr
    ELSE
      WRITE(6,'(1x,a,I2,a,/)') 'ERROR: gribtime in wrong length, ',itmp,  &
                               'It should be 13 or 14.'
      iret = -1
      RETURN
    END IF

    myr = MOD( iyr, 100 )
    IF (myr == 0) myr = 100

    INQUIRE(FILE=gribfile(1:grbflen), EXIST = fexist )

    IF( fexist ) THEN
      gribfile_new = gribfile
      grbflen_new = grbflen
      GO TO 200
    ENDIF

    INQUIRE(FILE=trim(gribfile(1:grbflen))//'.Z', EXIST = fexist )
    IF( fexist ) THEN
      CALL unixcmd( 'cp '//trim(gribfile(1:grbflen))//'.Z tem_ext_file.Z' )
      INQUIRE(FILE='tem_ext_file', EXIST = fexist )
      IF( fexist ) call unixcmd( '/bin/rm tem_ext_file' )
      CALL uncmprs( 'tem_ext_file.Z' )
      gribfile_new = 'tem_ext_file'
      grbflen_new = 12
      GO TO 200
    END IF

    INQUIRE(FILE=trim(gribfile(1:grbflen))//'.gz', EXIST = fexist )
    IF( fexist ) THEN
      CALL unixcmd( 'cp '//trim(gribfile(1:grbflen))//'.gz tem_ext_file.gz' )
      INQUIRE(FILE='tem_ext_file', EXIST = fexist )
      IF( fexist ) call unixcmd( '/bin/rm tem_ext_file' )
      CALL uncmprs( 'tem_ext_file.gz' )
      gribfile_new = 'tem_ext_file'
      grbflen_new = 12
      GO TO 200
    END IF

    WRITE(6,'(/1x,a,/1x,a/)')
    WRITE(6,'(1x,a,2(/,1x,a),/)') 'GRIB file '//gribfile(1:grbflen)//     &
        ' or its compressed version not found.',                          &
        'Please check option extdfmt and its instruction.',               &
        'Program warning in RDNMCGRB.'
    IRET = 1

! Note other errors go to 999 instead of 998.  We can't go to 999, as
! grbunit hasn't been allocated yet.

    GOTO 998

    200 CONTINUE

    WRITE(6,'(8A)') cmd, 'Opening GRIB ', gribfile(1:grbflen)

    opnmod = 'r'

    CALL gsize ( stdio_size )
    ALLOCATE(grbunit(stdio_size), STAT = istatus)
    CALL check_alloc_status(istatus, "grbunit")

    CALL gopen ( grbunit, gribfile_new, grbflen_new, opnmod, iret )

    IF  ( iret /= 0 )  THEN
    WRITE(6,'(a,a,/,a,i4)')                                             &
        'Error in opening file ',gribfile_new,                          &
        'C function gopen return status: ',iret
!     DEALLOCATE(grbunit)
!     RETURN
      IRET = 1
      GOTO 999
    END IF
!
!-----------------------------------------------------------------------
!
!  Get the size of the GRIB file
!
!-----------------------------------------------------------------------
!
    fsize = 0
    CALL gseek( grbunit, fsize, 2, iret )
    WRITE (6,'(2A,I10,A)') cmd, 'Size of GRIB file is', fsize,' bytes'
!
!-----------------------------------------------------------------------
!
!  Scan the GRIB file to count records
!
!-----------------------------------------------------------------------
!
    nrecs = 0
    CALL grbscan( grbunit, fsize, nprods, rcstart, rcbytes, nrecs )
    IF ( nrecs > nprods ) THEN
      WRITE (6,'(a/a,i6,a,i6/a/a)')                                       &
          'The actual number of products exceeded the expected number.',  &
          'Number of assumed products: ', nprods,                         &
          'Number of actual products: ', nrecs,                           &
          'Program stopped in RDNMCGRB. Please change NPRODS in ',        &
          'file gribcst.inc and re-run the program.'
      CALL gclose( grbunit, iret )
      CALL arpsstop("GRIB file problem",1)
    ELSE IF ( nrecs <= 0 ) THEN
      WRITE (6,'(a,a,a)')                                                 &
          'No GRIB message found in file ',gribfile_new,                  &
          'Program stopped in RDNMCGRB.'
      CALL gclose( grbunit, iret )
      CALL arpsstop("GRIB file problem",1)
    END IF
!
!-----------------------------------------------------------------------
!
!  Read the GRIB file and fill in the variable arrays
!
!-----------------------------------------------------------------------
!
  IF( n2dvars > n2dvs .OR. n2dlvtps > n2dlvt ) THEN
    WRITE(6,'(a,i3,a,i3,a)') 'Array var_nr2d(',n2dvs,',',n2dlvt,        &
                             ') declared in gribcst.inc too small.'
    WRITE(6,'(a,i3,a,i3,a)') 'It needs to be at least ',n2dvars,'x',n2dlvtps,'.'
    WRITE(6,'(a)') 'Program stopped in RDNMCGRB.'
    CALL arpsstop("GRIB file problem",1)
  ENDIF

  DO n=1,n2dvars
    DO m=1,n2dlvtps
      var_nr2d(n,m) = 0
    END DO
  END DO

  IF( n3dvars > n3dvs .OR. n3dlvtps > n3dlvt ) THEN
    WRITE(6,'(a,i3,a,i3,a)') 'Array var_nr3d(',n3dvs,',',n3dlvt,        &
                             ') declared in gribcst.inc too small.'
    WRITE(6,'(a,i3,a,i3,a)') 'It needs to be at least ',n3dvars,'x',n3dlvtps,'.'
    WRITE(6,'(a)') 'Program stopped in RDNMCGRB.'
    CALL arpsstop("GRIB file problem",1)
  ENDIF

  DO m=1,n3dlvtps
    DO n=1,n3dvars
      var_nr3d(n,m) = 0
    END DO
  END DO

  grdflg = .false.

  DO l=1,nrecs

    DO i=1,ipdsz
      ipds(i) = 0
    END DO

    DO i=1,igdsz
      igds(i) = 0
    END DO

    CALL gseek( grbunit, rcstart(l), 0, iret )
    CALL gread( grbunit, mgrib, rcbytes(l), iret )

    IF ( iendn == 1 ) THEN     ! little endian machine
      nwrd = (rcbytes(l)+3)/4
      CALL swap4byte( mgrib, nwrd )
    END IF

    CALL w3fi63( mgrib, ipds, igds, ibms, rcdata, mptrs, iret )

!    print '(1x,2a,I0,a,I0,5I5,F15.3)',cmd, 'PDS# ',L,'/',nrecs,         &
!                               igds(2),ipds(3),ipds(5),ipds(6),ipds(7),rcdata(1)
    IF ( ipds(3) /= gridtyp ) THEN
      IF (verbose) WRITE (6,'(/A,I5,2(A,I4))')                          &
         'Ignoring GRIB record', l, ': Grid ID is ipds(3) =',ipds(3),   &
         ', expected gridtyp=',gridtyp
      CYCLE
    END IF

    mlvtp  = 0
    chklev = 0
    DO m=1,n3dlvtps
      IF ( ipds(6) == levtyp3d(m) ) THEN
        chklev = 1
        mlvtp  = m
        GO TO 30
      END IF
    END DO

    DO m=1,n2dlvtps
      IF ( ipds(6) == levtyp2d(m) ) THEN
        IF (levtyp2d(m) == 105) THEN  ! T_2m, QV_2m, U_10m, V_10m
          SELECT CASE (ipds(5))
          CASE (11, 51, 52)           ! Must at 2m
            IF (ipds(7) /= 2) EXIT
          CASE (33, 34)               ! Must at 10m
            IF (ipds(7) /= 10) EXIT
          END SELECT
        END IF
        chklev = 2
        mlvtp  = m
        EXIT
      END IF
    END DO

    30      CONTINUE

    IF ( chklev == 0 ) THEN
      IF (verbose) &
        WRITE (6,'(/A,I5,A,I4)') 'Ignoring GRIB record', l, &
              ': Type of level not expected, ipds(6) =',ipds(6)
      CYCLE
    END IF

    IF ( ipds(16) /= 0 .AND. ipds(16) /= 10 .AND.           &
         ipds(16) /= 1 .AND. ipds(16) /= 4 ) THEN
      IF (verbose) &
        WRITE(6,'(/A,I5,A,I5)') 'Ignoring GRIB record', l,              &
            ': average, or time difference - ',ipds(16)
      CYCLE
    ELSE IF ( ipds(8)  /= myr  .OR. ipds(9)  /= imo .OR.                &
              ipds(10) /= iday .OR. ipds(11) /= ihr .OR.                &
              ipds(14) /= fhr ) THEN
      IF (ipds(16) == 4 .AND. ipds(15) == fhr) THEN
        !  good here for accumulated rain
      ELSE
        IF (verbose) WRITE (6,'(/A,I5,A,5I2.2,"f",I3.3,3A)')            &
              'Ignoring GRIB record', l, ': time of record (',          &
              ipds(21)-1,(ipds(j),j=8,11),ipds(14),                     &
              ') does not match expected (', gribtime, ")"
        CYCLE
      END IF
    END IF

    gdsflg = IAND(ipds(4),128)/128
    IF ( .NOT. grdflg ) THEN
      CALL grbgrid(gridtyp, gdsflg, igds,                               &
                   gridesc, iproj_grb, gthin,                           &
                   ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,       &
                   pi_grb,pj_grb,ipole, di_grb,dj_grb,                  &
                   latsw,lonsw, latne,lonne,                            &
                   latrot,lonrot,angrot,                                &
                   latstr,lonstr,facstr,                                &
                   lattru1,lattru2,lontrue,                             &
                   scanmode, iscan,jscan,kscan,                         &
                   ires,iearth,icomp,                                   &
                   jpenta,kpenta,mpenta,ispect,icoeff,                  &
                   xp_grb,yp_grb, xo_grb,yo_grb,zo_grb)

      WRITE (6,'(2A)') cmd, gridesc

      grdflg = .true.
    END IF

    IF ( igds(1) /= mproj_grb ) THEN
      IF (verbose) WRITE(6,'(/A,I5,2(A,I4))')                           &
         'Ignoring GRIB record', l,': Map projection is igds(1) =', igds(1), &
         ', expected mproj_grb =', mproj_grb
      CYCLE
    END IF

    IF ( igds(2) /= nx_ext .OR. igds(3) /= ny_ext) THEN
      IF (verbose) WRITE (6,'(a/a,i4,2x,a,i4/a,i4,2x,a,i4/a,I4,a)')     &
          'Horizontal dimension was not consistent,',                   &
          ' igds(2) = ', igds(2), ' igds(3) = ', igds(3),               &
          ' nx_ext  = ', nx_ext  , ' ny_ext  = ', ny_ext,               &
          'record: ',l,' skipped.'
      mismatch = mismatch + 1   ! skip mismatched record
      CYCLE                     ! stop only if all records are mismatched.
      !CALL arpsstop("GRIB data problem",1)
    END IF

    IF ( iscan == 0 ) THEN
      ileft  = 1
      iright = nx_ext
      iincr  = 1
    ELSE
      ileft  = nx_ext
      iright = 1
      iincr  = -1
    END IF

    IF ( jscan /= 0 ) THEN
      jleft  = 1
      jright = ny_ext
      jincr  = 1
    ELSE
      jleft  = ny_ext
      jright = 1
      jincr  = -1
    END IF

    IF ( kscan == 0 ) THEN
      ki = 1
      kj = nx_ext
    ELSE IF ( kscan == 1 ) THEN
      ki = ny_ext
      kj = 1
    ELSE
      IF (verbose) &
        WRITE(6,'(/A,I5,2(A,I4))') 'Ignoring GRIB record', l, &
                    'Unknown scanning mode, kscan = ', kscan
      CYCLE
    END IF

    IF (verbose) &
      WRITE (6,'(/A,I5,3(A,I5),$)') 'GRIB record', l, &
      ': param =', ipds(5), ', level type =', ipds(6), ', level =', ipds(7)

    IF ( chklev == 1 ) THEN           ! 3-D variable
      DO n=1,n3dvars
        !DO m=1,n3dlvtps   ! we already knew which level type it is,
                           ! then why do this loop?

        IF ( ipds(5) == var_id3d(n,mlvtp) ) THEN
          IF (verbose) WRITE (6,'(A)') '  accept 3D.'
          var_nr3d(n,mlvtp) = var_nr3d(n,mlvtp) + 1

          IF (levtyp3d(mlvtp) == 111 .OR. levtyp3d(mlvtp) == 112) THEN
!
!-----------------------------------------------------------------------
!
!  Insert the level and data values into the var_lev3d and var_grb3d
!  in the order that they occur in the data file (e.g. below surface
!  coordinates).
!
!  This data must be sorted.  When input is CONDUIT (LDM), data isn't
!  guaranteed to be in the proper order in the file!
!
!-----------------------------------------------------------------------
!
!              ilev = var_nr3d(n,mlvtp)

              ilev= 0
              DO k=1,nz_ext
                IF ( ipds(7) == var_lev3d(k,n,mlvtp) ) THEN
                  WRITE(6,*) 'Duplicate report, non-pressure coordinates'
                  var_nr3d(n,mlvtp) = var_nr3d(n,mlvtp) - 1
                  GOTO 555
                ENDIF
                IF ( ipds(7) < var_lev3d(k,n,mlvtp) .OR.     &
                     var_lev3d(k,n,mlvtp) < 0.0) THEN
                  ilev = k
                  EXIT
                END IF
              END DO

              IF ( ilev == 0 ) THEN
                PRINT *, 'couldnt locate level ',ipds(7)
                PRINT *, 'for var ID#',var_id3d(n,mlvtp)
                PRINT *, ' n = ', n, ' m = ', mlvtp
                PRINT *, ' levtyp3d = ',levtyp3d(mlvtp)
!               STOP
                CALL arpsstop("GRIB data problem",1)
              END IF
            ELSE
!
!-----------------------------------------------------------------------
!
!  Insert the level and data values into the var_lev3d and var_grb3d
!  arrays in order of decreasing level value,
!
!  For example if the level type is pressure levels the levels
!  will be sorted in order of decreasing pressure
!
!  ilev is the level at which to insert the current data
!
!-----------------------------------------------------------------------
!
              ilev= 0
              DO k=1,nz_ext
                IF ( ipds(7) == var_lev3d(k,n,mlvtp) ) then
                  WRITE(6,*) 'Duplicate report, pressure coordinates'
                  var_nr3d(n,mlvtp) = var_nr3d(n,mlvtp) - 1
                  GOTO 555
                ENDIF
                IF ( ipds(7) > var_lev3d(k,n,mlvtp) ) THEN
                  ilev = k
                  EXIT
                END IF
              END DO

              IF ( ilev == 0 ) THEN
                PRINT *, 'couldnt locate level ',ipds(7)
                PRINT *, 'for var ID#',var_id3d(n,mlvtp)
                PRINT *, ' n = ', n, ' m = ', mlvtp
!               STOP
                CALL arpsstop("GRIB data problem",1)
              END IF

            ENDIF

            ! slide all the following levels down one to make room

            DO k=nz_ext,ilev+1, -1
              var_lev3d(k,n,mlvtp)= var_lev3d(k-1,n,mlvtp)
              DO j=jleft,jright,jincr
                DO i=ileft,iright,iincr
                  var_grb3d(i,j,k,n,mlvtp)= var_grb3d(i,j,k-1,n,mlvtp)
                END DO
              END DO
            END DO

!
!   now insert the data at the appropriate level
!
            var_lev3d(ilev,n,mlvtp) = ipds(7)
            DO j=jleft,jright,jincr
              DO i=ileft,iright,iincr
                ij = ki*iincr*(i-ileft) + kj*jincr*(j-jleft) + 1
                var_grb3d(i,j,ilev,n,mlvtp) = rcdata(ij)
              END DO
            END DO
          END IF

          555  CONTINUE
        !END DO
      END DO
    ELSE IF ( chklev == 2 ) THEN       ! 2-D variable
      IF ( ipds(16) == 0 .OR. ipds(16) == 10 .OR.               &
           ipds(16) == 1 .OR. ipds(16) == 4 ) THEN
                                       ! Don't use average,
                                       ! or time difference fields.
        DO n=1,n2dvars
!          DO m=1,n2dlvtps
            IF ( ipds(5) == var_id2d(n,mlvtp) ) THEN
              IF (verbose) WRITE (6,'(a)') '  accept 2D.'
              var_nr2d(n,mlvtp) = var_nr2d(n,mlvtp) + 1
              DO i=ileft,iright,iincr
                DO j=jleft,jright,jincr
                  ij = ki*iincr*(i-ileft) + kj*jincr*(j-jleft) + 1
                  var_grb2d(i,j,n,mlvtp) = rcdata(ij)
                END DO
              END DO
            END IF
!          END DO
        END DO
      END IF
    END IF

  END DO
  IF (mismatch >= nrecs) THEN
    WRITE(6,'(1x,a,/,1x,2(a,I5),/,1x,2(a,I5),a)')                                      &
      'There is no record in the GRIB file that matches the desired dimension size.',  &
      '  expected: nx_ext  = ',nx_ext, ', ny_ext  = ',ny_ext,                          &
      '  found:    igds(2) = ',igds(2),', igds(3) = ',igds(3),'.'
    GOTO 999
  END IF
!
!-----------------------------------------------------------------------
!
!  Set good status
!
!-----------------------------------------------------------------------
!
  iret = 0

  999 CONTINUE

  ! Don't trash out "iret".

  CALL gclose( grbunit, itmp )
  DEALLOCATE( grbunit )

  998 CONTINUE

  END IF         ! IF (myproc == 0)

  CALL mpupdatei(iproj_grb,1)
  CALL mpupdatei(gthin,1)
  CALL mpupdatei(ni_grb,1)
  CALL mpupdatei(nj_grb,1)
  CALL mpupdatei(np_grb,1)
  CALL mpupdatei(nk_grb,1)
  CALL mpupdater(zk_grb,nz_ext)
  CALL mpupdatei(npeq,1)
  CALL mpupdater(nit,nz_ext)
  CALL mpupdater(pi_grb,1)
  CALL mpupdater(pj_grb,1)
  CALL mpupdatei(ipole,1)
  CALL mpupdater(di_grb,1)
  CALL mpupdater(dj_grb,1)
  CALL mpupdater(latsw,1)
  CALL mpupdater(lonsw,1)
  CALL mpupdater(latne,1)
  CALL mpupdater(lonne,1)
  CALL mpupdater(latrot,1)
  CALL mpupdater(lonrot,1)
  CALL mpupdater(angrot,1)
  CALL mpupdater(latstr,1)
  CALL mpupdater(lonstr,1)
  CALL mpupdater(facstr,1)
  CALL mpupdater(lattru1,1)
  CALL mpupdater(lattru2,1)
  CALL mpupdater(lontrue,1)
  CALL mpupdatei(scanmode,1)
  CALL mpupdatei(iscan,1)
  CALL mpupdatei(jscan,1)
  CALL mpupdatei(kscan,1)
  CALL mpupdatei(ires,1)
  CALL mpupdatei(iearth,1)
  CALL mpupdatei(icomp,1)
  CALL mpupdatei(jpenta,1)
  CALL mpupdatei(kpenta,1)
  CALL mpupdatei(mpenta,1)
  CALL mpupdatei(ispect,1)
  CALL mpupdatei(icoeff,1)
  CALL mpupdater(xp_grb,1)
  CALL mpupdater(yp_grb,1)
  CALL mpupdater(xo_grb,1)
  CALL mpupdater(yo_grb,1)
  CALL mpupdater(zo_grb,1)
  CALL mpupdater(rcdata,nx_ext*ny_ext)
  CALL mpupdater(var_grb2d,nx_ext*ny_ext*n2dvs*n2dlvt)
  CALL mpupdater(var_grb3d,nx_ext*ny_ext*nz_ext*n3dvs*n3dlvt)
  CALL mpupdatei(var_lev3d,nz_ext*n3dvs*n3dlvt)
  CALL mpupdatei(iret,1)

  CALL mpupdatei(n2dvars,1)
  CALL mpupdatei(n2dlvtps,1)
  CALL mpupdatei(var_nr2d,n2dvs*n2dlvt)

  CALL mpupdatei(n3dvars,1)
  CALL mpupdatei(var_nr3d,n3dvs*n3dlvt)

  RETURN
END SUBROUTINE rdnmcgrb
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GRBSCAN                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grbscan(grbunit,fsize,nprods,rcstart,rcbytes,nrecs)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is to scan the NMC GRIB files.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/03/1996 Initialized.
!
!  MODIFICATIONS:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    grbunit  Unit of GRIB file opened to scan for GRIB records
!    nprods   Maxmum record number
!    fsize    Size of the GRIB file
!
!  OUTPUT:
!
!    rcstart  Starting position (bytes) for each record
!    rcbytes  Record length in bytes
!    nrecs    Record number actually scanned
!    iret     Return flag
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: grbunit
  INTEGER :: fsize
  INTEGER :: nprods

  INTEGER :: rcbytes(nprods)  ! record length in bytes
  INTEGER :: rcstart(nprods)  ! record starting byte in a GRIB file

  INTEGER :: nrecs
  INTEGER :: iret
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nbytes

  CHARACTER (LEN=1) :: buffer(8)
!
!-----------------------------------------------------------------------
!
!  External function
!
!-----------------------------------------------------------------------
!
  INTEGER :: char2i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nbytes = 0
  nrecs  = 0

  100   IF ( nbytes < fsize ) THEN
  CALL gseek( grbunit,nbytes,0,iret )
  IF ( iret /= 0 ) GO TO 200

  CALL gread( grbunit,buffer,8,iret )
  IF ( iret /= 8 ) GO TO 200

  IF ( char2i(buffer(1)) /= 71 .OR.     & ! test ASCII 'G'
       char2i(buffer(2)) /= 82 .OR.     & ! test ASCII 'R'
       char2i(buffer(3)) /= 73 .OR.     & ! test ASCII 'I'
       char2i(buffer(4)) /= 66 .OR.     & ! test ASCII 'B'
       char2i(buffer(8)) /= 1 )      THEN ! test Edition #, 1
    nbytes = nbytes + 1
    GO TO 100
  END IF

  nrecs = nrecs + 1
  IF ( nrecs <= nprods ) THEN
    rcstart(nrecs) = nbytes
    rcbytes(nrecs) = char2i(buffer(5))*65536                            &
                   + char2i(buffer(6))*256                              &
                   + char2i(buffer(7))
    nbytes = nbytes + rcbytes(nrecs)
  END IF
  !PRINT *, 'nrecs = ',nrecs, rcbytes(nrecs),nbytes

  GO TO 100

  END IF

  200   RETURN

END SUBROUTINE grbscan
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE rdgrbdims                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdgrbdims(gribfile,grbflen,ni_grb,nj_grb,iret)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is to read dimension size from NMC GRIB #1 file.
!
!  The decoder of GRIB is from NMC, ftp:nic.fb4.noaa.gov/pub
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/12/2005 Initialized.
!
!  MODIFICATIONS:
!  08/17/2007 Yunheng Wang
!  Interface changed.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    gribfile      GRIB file name
!    grbflen       Length of GRIB file name
!
!  OUTPUT:
!
!    ni_grb        Grid size in west-east direction
!    nj_grb        Grid size in south-north direction
!
!    iret          Return flag
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: gribfile
  INTEGER,          INTENT(IN)  :: grbflen

  INTEGER,          INTENT(OUT) :: ni_grb       ! Number of points along x-axis
  INTEGER,          INTENT(OUT) :: nj_grb       ! Number of points along y-axis

  INTEGER,          INTENT(OUT) :: iret         ! Return flag
!
!-----------------------------------------------------------------------
!
!  Temporary arrays to read GRIB file
!
!-----------------------------------------------------------------------
!
  INTEGER :: nrecs        ! number of records

  INTEGER :: iproj        ! Map projection indicator

  INTEGER :: fsize        ! Size of GRIB file
  INTEGER :: stdio_size   ! Size of "struct file" (C code in gribio_c.c).

  CHARACTER(LEN=1), ALLOCATABLE :: grbunit(:) ! C structure, size stdio_size
  LOGICAL (KIND=1), ALLOCATABLE :: ibms(:)    ! BMS logical array for data bit map
  REAL,             ALLOCATABLE :: rcdata(:)  ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Misc. Local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=1)   :: opnmod

  INTEGER             :: iendn,itypec,wdlen, nwrd
  INTEGER             :: l, i

  LOGICAL,           PARAMETER :: verbose = .FALSE.
  CHARACTER(LEN=11), PARAMETER :: cmd     = 'rdgrbdims: '

  INCLUDE 'gribcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(ibms(nbufsz/4),   STAT = iret)
  CALL check_alloc_status(iret, 'rdgrbdims:ibms')
  ALLOCATE(rcdata(nbufsz/4), STAT = iret)
  CALL check_alloc_status(iret, 'rdgrbdims:rcdata')
!
!-----------------------------------------------------------------------
!
!  Find the type of endian, type of character set and length of
!  machine word
!
!  iendn     -  Integer for big-endian or little-endian
!               = 0   big-endian
!               = 1   little-endian
!               = 2   cannot compute
!  itypec    -  Integer for type of character set
!               = 0   ascii  character set
!               = 1   ebcdic character set
!               = 2   not ascii or ebcdic
!  wdlen     -  Integer for words size of computer in bytes
!               = 4   for 32 bit (4 bytes) computers
!               = 8   for 64 bit (8 bytes) computers
!
!-----------------------------------------------------------------------
!
  CALL findvartype(iendn,itypec,wdlen)
  IF ( iendn == 2 ) THEN
    WRITE(6,'(a)') 'Unknown endian type. Stopped in rdnmcgrb'
    iret = -1
    GOTO 999
  END IF

  IF ( itypec == 2 ) THEN
    WRITE(6,'(a)') 'Unknown character set. Stopped in rdnmcgrb'
    iret = -1
    GOTO 999
  END IF
!
!-----------------------------------------------------------------------
!
!  Open the GRIB file
!
!-----------------------------------------------------------------------
!
    WRITE(6,'(8A)') cmd, 'Opening GRIB ', TRIM(gribfile)

    opnmod = 'r'

    CALL gsize ( stdio_size )
    ALLOCATE(grbunit(stdio_size), STAT = iret)
    CALL check_alloc_status(iret, "rdgrbdim:grbunit")

    CALL gopen ( grbunit, gribfile, grbflen, opnmod, iret )

    IF  ( iret /= 0 )  THEN
      WRITE(6,'(a,a,/,a,i4)') 'Error in opening file ',gribfile,    &
                              'C function gopen return status: ',iret
      DEALLOCATE(grbunit)
      GOTO 999
    END IF
!
!-----------------------------------------------------------------------
!
!  Get the size of the GRIB file
!
!-----------------------------------------------------------------------
!
    fsize = 0
    CALL gseek( grbunit, fsize, 2, iret )
    WRITE (6,'(2A,I10,A)') cmd, 'Size of GRIB file is', fsize,' bytes'
!
!-----------------------------------------------------------------------
!
!  Scan the GRIB file to count records
!
!-----------------------------------------------------------------------
!
    nrecs = 0
    CALL grbscan( grbunit, fsize, nprods, rcstart, rcbytes, nrecs )
    IF ( nrecs > nprods ) THEN
      WRITE (6,'(a/a,i6,a,i6/a/a)')                                     &
        'The actual number of products exceeded the expected number.',  &
        'Number of assumed products: ', nprods,                         &
        'Number of actual products: ', nrecs,                           &
        'Program stopped in RDNMCGRB. Please change NPRODS in ',        &
        'file gribcst.inc and re-run the program.'
      CALL gclose( grbunit, iret )
      iret = -3
      GOTO 999
    ELSE IF ( nrecs <= 0 ) THEN
      WRITE (6,'(a,a,a)')                                               &
        'No GRIB message found in file ',gribfile,                      &
        'Program stopped in RDNMCGRB.'
      CALL gclose( grbunit, iret )
      iret = -3
      GOTO 999
    END IF
!
!-----------------------------------------------------------------------
!
!  Read the first record in GRIB file
!
!-----------------------------------------------------------------------
!
    DO l = 1,nrecs

      DO i=1,ipdsz
        ipds(i) = 0
      END DO
      DO i=1,igdsz
        igds(i) = 0
      END DO

      CALL gseek( grbunit, rcstart(l), 0, iret )
      CALL gread( grbunit, mgrib, rcbytes(l), iret )

      IF ( iendn == 1 ) THEN       ! little endian machine
        nwrd = (rcbytes(l)+3)/4
        CALL swap4byte( mgrib, nwrd )
      END IF

      CALL w3fi63( mgrib, ipds, igds, ibms, rcdata, mptrs, iret )

      iproj = IABS(igds(1))
      IF (iproj ==  0 .OR. iproj == 3) THEN
        ni_grb = igds(2)
        nj_grb = igds(3)
      ELSE
        WRITE(6,'()') 'Unknow map projection: ',iproj
        iret = -5
        EXIT
      END IF

      iret = 0  ! set good status
      EXIT

    END DO

!
!-----------------------------------------------------------------------
!
!  close grib file
!
!-----------------------------------------------------------------------
!
    CALL gclose( grbunit, iret )
    DEALLOCATE( grbunit )

  999 CONTINUE
  DEALLOCATE(ibms,rcdata)

  RETURN
END SUBROUTINE rdgrbdims
