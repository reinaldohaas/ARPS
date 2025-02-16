!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDNMCGRB2                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE rdnmcgrb2(gribfile,grbflen, gribtime,                        &
           gridesc, iproj_grb, gthin,                                   &
           ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,               &
           pi_grb,pj_grb,ipole, di_grb,dj_grb,                          &
           latsw,lonsw, latne,lonne,                                    &
           latrot,lonrot,angrot,                                        &
           latstr,lonstr,facstr,                                        &
           lattru1,lattru2,lontrue,                                     &
           scanmode, iscan,jscan,kscan,                                 &
           ires,iearth,icomp,                                           &
           jpenta,kpenta,mpenta,ispect,icoeff,                          &
           xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                         &
           nxgrb,nygrb,nzgrb,var_lev3d,var_grb3d,                       &
           var_grb2d,ibms,rcdata, iret)

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
!  November 1999 (Eric Kemp)
!  Added several arguments to subroutine call.
!
!  9 January 2001 (Eric Kemp)
!  Added proper Y2K bug fix conforming with ARPS 5.0.
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

  USE gribcst2

  IMPLICIT NONE

  INTEGER :: nxgrb,nygrb,nzgrb
  INTEGER :: var_lev3d(nzgrb,n3dvs,n3dlvt)
     ! Levels (hybrid) for each 3-D variable
  REAL :: var_grb3d(nxgrb,nygrb,nzgrb,n3dvs,n3dlvt) ! GRIB 3-D variables
  REAL :: var_grb2d(nxgrb,nygrb,n2dvs,n2dlvt) ! GRIB variables
  LOGICAL :: ibms(nxgrb*nygrb)      ! BMS logical array for data bit map
  REAL :: rcdata(nxgrb*nygrb)      ! temporary data array

  INTEGER :: grbflen
  CHARACTER (LEN=80) :: gribfile
  CHARACTER (LEN=13) :: gribtime
  CHARACTER (LEN=1) :: opnmod

!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: gdsflg       ! GDS indicator

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nzgrb)   ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nzgrb)   ! Number of x-points for thinned grid

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

!-----------------------------------------------------------------------
!
!  Temporary arrays to read GRIB file
!
!-----------------------------------------------------------------------

  INTEGER :: nrecs        ! number of records

  INTEGER :: fsize        ! Size of GRIB file
  INTEGER :: grbunit      ! I/O unit of GRIB file

!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,l,m,n, ij, ki, kj, ilev
  INTEGER :: iyr,imo,iday,ihr,fhr, myr

  INTEGER :: ileft, iright, iincr
  INTEGER :: jleft, jright, jincr

  INTEGER :: chklev, mlvtp
  INTEGER :: iendn,itypec,wdlen, nwrd

  LOGICAL :: grdflg           ! flag for grid information
  LOGICAL :: fexist

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

  CALL w3fi04(iendn,itypec,wdlen)
  IF ( iendn == 2 ) THEN
    WRITE(6,'(a)') 'Unknown endian type. Stopped in rdnmcgrb2'
    STOP
  END IF

  IF ( itypec == 2 ) THEN
    WRITE(6,'(a)') 'Unknown character set. Stopped in rdnmcgrb2'
    STOP
  END IF

!-----------------------------------------------------------------------
!
!   initialize the var_lev3d array so we can sort the records by
!   level values as we read them. We can't assume monotonic ordering
!   of the levels
!
!-----------------------------------------------------------------------

  DO k=1, nzgrb
    DO n=1, n3dvs
      DO m=1, n3dlvt
        var_lev3d(k,n,m)= -999999
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Open the GRIB file
!
!-----------------------------------------------------------------------

  READ (gribtime,'(i4,3i2,1x,i2)') iyr,imo,iday,ihr,fhr

!  myr = iyr - 1900
!arpscntl Y2K RLC per GMB/Donghai
! myr = MOD( iyr, 100 )

!EMK...Proper Y2K bug fix found in ARPS 5.0 package
  myr = MOD( iyr, 100 )   
  IF (myr == 0) THEN
    myr = 100
  END IF


  INQUIRE(FILE=gribfile(1:grbflen), EXIST = fexist )
  IF ( .NOT.fexist ) THEN
    WRITE(6,'(a)') 'GRIB file '//gribfile(1:grbflen)//                  &
        ' not found. Program stopped in RDNMCGRB2.'
    STOP
  END IF

  PRINT *, ' Opening GRIB file, gribfile= ',gribfile
  opnmod = 'r'
  CALL gopen  ( grbunit, gribfile, grbflen, opnmod, iret )
  IF  ( iret /= 0 )  THEN
    WRITE(6,'(a,a,/,a,i4)')                                             &
        'Error in opening file ',gribfile,                              &
        'C function gopen return status: ',iret
    RETURN
  END IF

!-----------------------------------------------------------------------
!
!  Get the size of the GRIB file
!
!-----------------------------------------------------------------------

  fsize = 0
  CALL gseek( grbunit, fsize, 2, iret )
  WRITE (6,'(a,i16,a)') 'The size of GRIB file is ',fsize,' bytes'

!-----------------------------------------------------------------------
!
!  Scan the GRIB file to count records
!
!-----------------------------------------------------------------------

  nrecs = 0
  CALL grbscan( grbunit, fsize, nprods, rcstart, rcbytes, nrecs )
  IF ( nrecs > nprods ) THEN
    WRITE (6,'(a/a,i6,a,i6/a/a)')                                       &
        'The actual number of products exceeded the expected number.',  &
        'Number of assumed products: ', nprods,                         &
        'Number of actual products: ', nrecs,                           &
        'Program stopped in RDNMCGRB2. Please change NPRODS in ',       &
        'file gribcst.inc and re-run the program.'
    CALL gclose( grbunit, iret )
    STOP
  ELSE IF ( nrecs <= 0 ) THEN
    WRITE (6,'(a,a,a)')                                                 &
        'No GRIB message found in file ',gribfile,                      &
        'Program stopped in RDNMCGRB2.'
    CALL gclose( grbunit, iret )
    STOP
  END IF

!-----------------------------------------------------------------------
!
!  Read the GRIB file and fill in the variable arrays
!
!-----------------------------------------------------------------------

  DO n=1,n2dvars
    DO m=1,n2dlvtps
      var_nr2d(n,m) = 0
    END DO
  END DO

  DO m=1,n3dlvtps
    DO n=1,n3dvars
      var_nr3d(n,m) = 0
    END DO
  END DO

  grdflg = .false.

  WRITE(6,*)'Unpacking ',nrecs,' records...'
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

    IF ( ipds(3) /= gridtyp ) THEN
!      write (6,'(a,i6,a/a,i6,a,i6)') 'Grid id of record, ',L,
!    :    'was inconsistant with expected,',
!    :    'ipds(3)=',ipds(3), ' gridtyp=',gridtyp
      CYCLE
    END IF

    DO m=1,n3dlvtps
      IF ( ipds(6) == levtyp3d(m) ) THEN
        chklev = 1
        mlvtp  = m
        GO TO 30
      ELSE
        chklev = 0
      END IF
    END DO

    DO m=1,n2dlvtps
      IF ( ipds(6) == levtyp2d(m) ) THEN
        chklev = 2
        mlvtp  = m
        EXIT
      ELSE
        chklev = 0
      END IF
    END DO

    30      CONTINUE

    IF ( chklev == 0 ) THEN
!      write (6,'(a,i4,a/a,i4)')
!    :    'Type of level of record, ',L,
!    :    ' was inconsistant with expected,',
!    :    ' ipds(6) =',ipds(6)
      CYCLE
    END IF

    IF ( ipds(8) /= myr .OR.ipds(9) /= imo.OR.                          &
           ipds(10) /= iday.OR.ipds(11) /= ihr.OR.                      &
           ipds(14) /= fhr ) THEN
!arpscntl RLC format
      WRITE (6,'(a,i6,a/a,I3.2,3i2.2,i2.2/2a)')                         &
          'Time of record, ',l, 'was inconsistant with expected,',      &
          ' pds time = ',(ipds(j),j=8,11),ipds(14),                     &
          ' gribtime = ', TRIM(gribtime)
      WRITE (6,'(A,I3.2,3i2.2,i2.2)') "Gribtime: ", myr, imo, iday, ihr, fhr
      CYCLE
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

      WRITE (6,'(a)') gridesc

      grdflg = .true.
    END IF

    IF ( igds(1) /= mproj_grb ) THEN
!      write (6,'(a,i3)')
!    :    'Map projection was not consistant, igds(1) = ', igds(1)
      CYCLE
    END IF

    IF ( igds(2) /= nxgrb .OR. igds(3) /= nygrb ) THEN
      WRITE (6,'(a/a,i4,2x,a,i4/a,i4,2x,a,i4/a)')                       &
          'Horizontal dimension was not consistant,',                   &
          ' igds(2) = ', igds(2), ' igds(3) = ', igds(3),               &
          ' nxgrb   = ', nxgrb  , ' nygrb   = ', nygrb,                 &
          'Program stopped in RDNMCGRB2.'
      STOP
    END IF

    IF ( iscan == 0 ) THEN
      ileft  = 1
      iright = nxgrb
      iincr  = 1
    ELSE
      ileft  = nxgrb
      iright = 1
      iincr  = -1
    END IF

    IF ( jscan /= 0 ) THEN
      jleft  = 1
      jright = nygrb
      jincr  = 1
    ELSE
      jleft  = nygrb
      jright = 1
      jincr  = -1
    END IF

    IF ( kscan == 0 ) THEN
      ki = 1
      kj = nxgrb
    ELSE IF ( kscan == 1 ) THEN
      ki = nygrb
      kj = 1
    ELSE
!      write (6,'(a,i3,a,i6)')
!    :    'Unknown scanning mode, kScan = ', kScan
      CYCLE
    END IF

    IF ( chklev == 1 ) THEN           ! 3-D variable
      DO n=1,n3dvars
        DO m=1,n3dlvtps

          IF ( ipds(5) == var_id3d(n,m) ) THEN
            var_nr3d(n,m) = var_nr3d(n,m) + 1

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

            ilev= 0
            DO k=1,nzgrb
              IF ( ipds(7) > var_lev3d(k,n,m) ) THEN
                ilev = k
                EXIT
              END IF
            END DO

            120           CONTINUE

            IF ( ilev == 0 ) THEN
              PRINT*,'couldnt locate level ',ipds(7)
              PRINT*,'for var ID#',var_id3d(n,m)
              STOP
            END IF

!   slide all the following levels down one to make room

            DO k=nzgrb, ilev+1, -1
              var_lev3d(k,n,m)= var_lev3d(k-1,n,m)
              DO j=jleft,jright,jincr
                DO i=ileft,iright,iincr
                  ij = ki*iincr*(i-ileft) + kj*jincr*(j-jleft) + 1
                  var_grb3d(i,j,k,n,m)= var_grb3d(i,j,k-1,n,m)
                END DO
              END DO
            END DO

!   now insert the data at the appropriate level

            var_lev3d(ilev,n,m) = ipds(7)
            DO j=jleft,jright,jincr
              DO i=ileft,iright,iincr
                ij = ki*iincr*(i-ileft) + kj*jincr*(j-jleft) + 1
                var_grb3d(i,j,ilev,n,m) = rcdata(ij)
              END DO
            END DO
          END IF
        END DO
      END DO
    ELSE IF ( chklev == 2 ) THEN       ! 2-D variable
      DO n=1,n2dvars
        DO m=1,n2dlvtps
          IF ( ipds(5) == var_id2d(n,m) ) THEN
            var_nr2d(n,m) = var_nr2d(n,m) + 1
            DO i=ileft,iright,iincr
              DO j=jleft,jright,jincr
                ij = ki*iincr*(i-ileft) + kj*jincr*(j-jleft) + 1
                var_grb2d(i,j,n,m) = rcdata(ij)
              END DO
            END DO
          END IF
        END DO
      END DO
    END IF

  END DO

!-----------------------------------------------------------------------
!
!  Set good status
!
!-----------------------------------------------------------------------

  iret = 0

  CALL gclose( grbunit, iret )

  RETURN
END SUBROUTINE rdnmcgrb2
