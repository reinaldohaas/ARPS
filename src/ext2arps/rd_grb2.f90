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
!

SUBROUTINE rdnmcgrb2(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,   &
           ix1,ix2,jy1,jy2,                                             &
           iproj_grb,nx_grb,ny_grb,dx_grb,dy_grb,                       &
           latsw,lonsw,lattru1,lattru2,lontrue,uvearth,                 &
           n2dvs, n3dvs, maxvar, nzsoilin_ext,                          &
           varids, var2dindx, var3dindx, var2dlvl, var3dlvl, var3dsoil, &
           var_grb2d, var_grb3d, lvldbg,                                &
           iret)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is to read and select given variables with GRIB
!  grid ID and level type from NMC GRIB files.
!
!  The decoder of GRIB2 is from NCEP.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/01/2006
!
!  MODIFICATIONS:
!  Keith Brewster (6/30/2011)
!  Changed logic for grid type 104 (sigma levels) which can be either
!  a 2D or 3D field.   Logic depends on var3dindx
!   var3dindx(m) > 0  3D Field
!   var3dindx(m) = 0  and var2dindx(m) > 0, 2D Field
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
  USE GRIB_MOD

  IMPLICIT NONE

  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN) :: nx_ext,ny_ext,nz_ext

  INTEGER, INTENT(IN) :: grbflen

  CHARACTER(LEN=*), INTENT(IN) :: gribfile
  CHARACTER(LEN=*), INTENT(IN) :: gribtime

  INTEGER, INTENT(IN)  :: ix1, ix2, jy1, jy2

  INTEGER, INTENT(OUT) :: iproj_grb    ! Map projection indicator
                                       ! Already converted to ARPS map definitions
  INTEGER, INTENT(OUT) :: nx_grb       ! Number of points along x-axis
  INTEGER, INTENT(OUT) :: ny_grb       ! Number of points along y-axis
  REAL,    INTENT(OUT) :: dx_grb       ! x-direction increment or grid length
  REAL,    INTENT(OUT) :: dy_grb       ! y-direction increment or grid length
  REAL,    INTENT(OUT) :: latsw        ! Latitude  of South West corner point
  REAL,    INTENT(OUT) :: lonsw        ! Longitude of South West corner point
  REAL,    INTENT(OUT) :: lattru1      ! Latitude (1st) at which projection is true
  REAL,    INTENT(OUT) :: lattru2      ! Latitude (2nd) at which projection is true
  REAL,    INTENT(OUT) :: lontrue      ! Longitude      at which projection is true
  INTEGER, INTENT(OUT) :: uvearth      ! = 0, Resolved u and v components of vector
                               ! quantities relative to easterly and northerly directions
                                       ! = 1, Resolved u and v components of vector quantities relative to the defined grid in the direction of increasing x and y (or i and j) coordinates, respectively

  INTEGER, INTENT(IN)  :: n2dvs, n3dvs, maxvar
  INTEGER, INTENT(IN)  :: nzsoilin_ext
  INTEGER, INTENT(IN)  :: varids(4,maxvar)
  INTEGER, INTENT(IN)  :: var2dindx(maxvar), var3dindx(maxvar)
  REAL,    INTENT(IN)  :: var2dlvl(n2dvs), var3dlvl(nz_ext), var3dsoil(nzsoilin_ext)
  REAL,    INTENT(OUT) :: var_grb2d(nx_ext,ny_ext,n2dvs)
  REAL,    INTENT(OUT) :: var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs)

  INTEGER, INTENT(IN)  :: lvldbg

  INTEGER, INTENT(OUT) :: iret         ! Return flag
!
!-----------------------------------------------------------------------
!
!  Temporary arrays to read GRIB file
!
!-----------------------------------------------------------------------
!
  INTEGER :: grbunit

  CHARACTER(LEN=1), ALLOCATABLE :: cgrib(:)

  INTEGER, PARAMETER :: maxlen = 32000

  TYPE(gribfield) :: gfld

!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,m,n,ij
  INTEGER :: ii, jj, ioff, joff
  INTEGER :: iyr,  imo,  iday,  ihr,  fhr

  LOGICAL             :: fexist
  INTEGER             :: grbflen_new
  CHARACTER (LEN=256) :: gribfile_new

  LOGICAL             :: verbose = .FALSE.
  CHARACTER(LEN=12), PARAMETER :: cmd     = 'rdnmcgrib2: '

  INTEGER :: iseek, icount, itot
  INTEGER :: lgrib, lgribin, lskip
  INTEGER :: currlen
  INTEGER :: listsec0(3), listsec1(13)
  INTEGER :: numfields, numlocal, maxlocal

  REAL    :: levelin

  INTEGER :: iscan, jscan, ijscan
  INTEGER :: ibgn, iinc, iend, jbgn, jinc, jend
  INTEGER :: imin, jmin

  INTEGER :: itemp

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0) THEN

    IF (lvldbg > 90) verbose = .TRUE.

!-----------------------------------------------------------------------
!
! Date and time expected
!
!-----------------------------------------------------------------------

    itemp = LEN_TRIM(gribtime)
    IF (itemp == 13) THEN
      READ (gribtime,'(i4,3i2,1x,i2)') iyr,imo,iday,ihr,fhr
    ELSE IF (itemp == 14 ) THEN
      READ (gribtime,'(i4,3i2,1x,i3)') iyr,imo,iday,ihr,fhr
    ELSE
      WRITE(6,'(1x,a,I2,a,/)') 'ERROR: gribtime in wrong length, ',itemp, &
                               'It should be 13 or 14.'
      iret = -1
      GOTO 9999
    END IF
!
!-----------------------------------------------------------------------
!
!  Open the GRIB file
!
!-----------------------------------------------------------------------
!
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
    WRITE(6,'(a)') 'GRIB file '//gribfile(1:grbflen)//                      &
        ' or its compressed version not found. Program stopped in RDNMCGRB2.'
    iret = -5
    GOTO 9999

    200 CONTINUE

    WRITE(6,'(1x,8A)') cmd, 'Opening GRIB ', gribfile(1:grbflen)

    CALL getunit( grbunit )

    CALL baopenr(grbunit,gribfile_new(1:grbflen_new),iret)
    IF  ( iret /= 0 )  THEN
      WRITE(6,'(a,a,/,a,i4)')                                             &
          'ERROR in opening file ',gribfile_new(1:grbflen_new),           &
          'function baopenr return status: ',iret
      GOTO 9999
    END IF

!-----------------------------------------------------------------------
!
! Loop over all records
!
!-----------------------------------------------------------------------

    itot    = 0     ! total field counter
    icount  = 0     ! total GRIB2 message counter
    iseek   = 0
    currlen = 0

  !  DO n = 1,maxvar
  !    WRITE(0,*) varids(:,n)
  !  END DO
  !  STOP

    iret = 0
    DO WHILE  (iret == 0)

      CALL skgb(grbunit,iseek,maxlen,lskip,lgrib)
      IF (lgrib == 0) EXIT

      IF (lgrib > currlen) THEN      ! Extend buffer size as needed
        IF (ALLOCATED(cgrib)) DEALLOCATE(cgrib)
        ALLOCATE(cgrib(lgrib), STAT = iret)
        currlen = lgrib
      END IF

      CALL baread(grbunit,lskip,lgrib,lgribin,cgrib)   ! read a GRIB message
      IF (lgrib /= lgribin) THEN
        WRITE(6,'(2(a,I6))')                                              &
           'ERROR: GRIB2 message size in error, required size = ',lgrib,  &
           ', read-in size = ',lgribin
        iret = -3
        GOTO 9999
!        CALL arpsstop('ERROR in rdnmcgrb2 when calling baread.',1)
      END IF
      icount = icount + 1
      iseek = lskip+lgrib            ! next GRIB message starting location

      CALL gb_info(cgrib,lgrib,listsec0,listsec1,                         &
                   numfields,numlocal,maxlocal,iret)
      IF (iret /= 0) THEN
        WRITE(6,'(a,I3)') 'ERROR: quering GRIB2 message = ',iret
        GOTO 9999
!        CALL arpsstop('ERROR in rdnmcgrb2 when calling gb_info.',1)
      END IF
      itot = itot+numfields

      IF (icount == 1) THEN          ! Once per file for date and map projection value
        CALL gf_getfld(cgrib,lgrib,1,.TRUE.,.TRUE.,gfld,iret)
        IF (iret /= 0) THEN
          WRITE(6,'(a,I3)') 'ERROR: Decoding GRIB2 message = ',iret
          GOTO 9999
!          CALL arpsstop('ERROR in rdnmcgrb2 when calling gf_getfld.',1)
        END IF

        IF (verbose)  &
          WRITE(6,'(1x,a,I4.4,4(a,I2.2),/)') 'GRIB date-time found is: ', &
              gfld%idsect(6),'-',gfld%idsect(7),'-',gfld%idsect(8),'_',   &
              gfld%idsect(9),'f',gfld%ipdtmpl(9)

        IF (iyr  /= gfld%idsect(6) .OR. imo /= gfld%idsect(7) .OR.      &
            iday /= gfld%idsect(8) .OR. ihr /= gfld%idsect(9) .OR.      &
            ( gfld%ipdtnum == 0 .AND. fhr  /= gfld%ipdtmpl(9) ) ) THEN
          WRITE(6,'(a,/,7x,2a,/,7x,a,I4.4,3I2.2,a,I2.2)')               &
              'ERROR: GRIB data is in wrong date and time',             &
              'Expected: ',gribtime, '   Found: ',gfld%idsect(6),       &
              gfld%idsect(7),gfld%idsect(8),gfld%idsect(9),             &
              'f',gfld%ipdtmpl(9)
          iret = -4
          GOTO 9999
!          CALL arpsstop('Wrong grib file?',1)
        END IF

        itemp = gfld%ipdtmpl(5)         ! indicator of model
        IF (lvldbg > 50) THEN
          IF (gfld%idsect(1) == 7) THEN
            IF (itemp == 83 .OR. itemp == 84) THEN
              WRITE(6,'(1x,a,/)') 'NCEP NAM model'
            ELSE IF (itemp == 81 .OR. itemp == 96) THEN
              WRITE(6,'(1x,a,/)') 'NCEP GFS model'
            ELSE
              WRITE(6,'(1x,a,I3,/)') 'Unkown model from NCEP, model = ',itemp
            END IF
          ELSE
            WRITE(6,'(1x,a,I2,a,I2,/)')                                   &
              'Unknown model and orig center, center = ',gfld%idsect(1),  &
              ' model = ',itemp
          END IF
        END IF

        IF (gfld%igdtnum == 0) THEN       ! Lat/Lon grid
          iproj_grb = 4
          nx_grb    = gfld%igdtmpl(8)
          ny_grb    = gfld%igdtmpl(9)
          dx_grb    = gfld%igdtmpl(17)/1.0E6  ! in degree
          dy_grb    = gfld%igdtmpl(18)/1.0E6  ! in degree
          iscan     = IAND(gfld%igdtmpl(19),128)/128  ! get the first bit
                      ! = 0, +x(+i) direction, = 1, -x(-i) direction
          jscan     = IAND(gfld%igdtmpl(19),64) / 64   ! get the second bit
                      ! = 0, -y(-j) direction, = 1, +y(+j) direction
          ijscan    = IAND(gfld%igdtmpl(19),32) / 32   ! get the third bit
                      ! = 0, x direction first, = 1, y direction first

          uvearth   = IAND(gfld%igdtmpl(14),8)  /  8   ! get the fifth bit, big-endian
          IF (jscan == 1) THEN          ! +y(+j) direction
            latsw     = gfld%igdtmpl(12)/1.0E6
            lonsw     = gfld%igdtmpl(13)/1.0E6
          ELSE                          ! -y(-j) direction
            latsw     = gfld%igdtmpl(15)/1.0E6
            lonsw     = gfld%igdtmpl(13)/1.0E6
          END IF
        ELSE IF (gfld%igdtnum == 30) THEN ! Lambert Conformal Grid
          iproj_grb   = 2
          nx_grb      = gfld%igdtmpl(8)
          ny_grb      = gfld%igdtmpl(9)
          lontrue     = gfld%igdtmpl(14)/1.0E6
          lattru1     = gfld%igdtmpl(19)/1.0E6
          lattru2     = gfld%igdtmpl(20)/1.0E6
          dx_grb      = gfld%igdtmpl(15)/1.0E3
          dy_grb      = gfld%igdtmpl(16)/1.0E3
          latsw       = gfld%igdtmpl(10)/1.0E6
          lonsw       = gfld%igdtmpl(11)/1.0E6
          iscan     = IAND(gfld%igdtmpl(18),128)/128  ! get the first bit
                      ! = 0, +x(+i) direction, = 1, -x(-i) direction
          jscan     = IAND(gfld%igdtmpl(18),64) / 64   ! get the second bit
                      ! = 0, -y(-j) direction, = 1, +y(+j) direction
          ijscan    = IAND(gfld%igdtmpl(18),32) / 32   ! get the third bit
                      ! = 0, x direction first, = 1, y direction first

        ELSE IF (gfld%igdtnum == 20) THEN ! Polar-Stereographic Grid
          iproj_grb   = 1
          nx_grb      = gfld%igdtmpl(8)
          ny_grb      = gfld%igdtmpl(9)
          lattru1     = 60.
          lattru2     = 91.
          latsw       = gfld%igdtmpl(10)
          lonsw       = gfld%igdtmpl(11)

        ELSE
           WRITE(6,'(1x,a,I3)') 'Unkown projection: ',gfld%igdtnum
           iret = -4
           GOTO 9999
!           CALL arpsstop('Unknown map project in GRIB2 file.',1)
        END IF

        !iinc = (-1)**iscan
        !ibgn = (nx_ext)**iscan
        !iend = ibgn + (nx_ext-1)*iinc
        !
        !jinc = (-1)*(-1)**jscan
        !jend = (ny_ext)**jscan
        !jbgn = jend - (ny_ext-1)*jinc
        !
        !IF (nx_ext /= nx_grb .OR. ny_ext /= ny_grb) THEN
        !  WRITE(6,'(1x,a,/,8x,a,I3,a,I3,/,8x,a,I3,a,I3)')                 &
        !     'ERROR: Dimension size inconsistent',                        &
        !     'Expected: nx_ext = ',nx_ext,', ny_ext = ',ny_ext,           &
        !     'Found   : nx_grb = ',nx_grb,', ny_grb = ',ny_grb
        !  iret = -4
        !  GOTO 9999
        !  !CALL arpsstop('Inconsistent dimension size in GRIB2 file.',1)
        !END IF

        IF (iscan == 0) THEN
          iinc = 1
          ibgn = MAX(1,ix1)
          iend = MIN(nx_grb,ix2)
          ioff = 0
          imin = ibgn-1
        ELSE
          iinc = -1
          ibgn = MIN(nx_grb,ix2)
          iend = MAX(1,ix1)
          ioff = nx_grb+1
          imin = iend-1
        END IF

        IF (jscan == 1) THEN
          jinc = 1
          jbgn = MAX(1,jy1)
          jend = MIN(ny_grb,jy2)
          joff = -1
          jmin = jbgn-1
        ELSE
          jinc = -1
          jbgn = MIN(ny_grb,jy2)
          jend = MIN(1,jy1)
          joff = ny_grb
          jmin = jend-1
        END IF

        IF (nx_ext /= ABS(iend-ibgn)+1 .OR. ny_ext /= ABS(jend-jbgn)+1 ) THEN
          WRITE(6,'(1x,a,/,8x,a,I3,a,I3,/,8x,a,I3,a,I3)')                 &
             'ERROR: Dimension size inconsistent',                        &
             'Expected: nx_ext = ',nx_ext,', ny_ext = ',ny_ext,           &
             'Found   : nx_grb = ',nx_grb,', ny_grb = ',ny_grb
          iret = -4
          GOTO 9999
        END IF

        IF (verbose) THEN
          WRITE(6,'(1x,2a)') 'M_No  StartBytes F_No. ',                   &
               'Field_No.  (Dis, Cat, Par) LevelType   Levelval  GRB_Array'
          WRITE(6,'(1x,2a)') '===== ========== ===== ',                   &
               '========== =============== ========== ========== =========='
        END IF

        CALL gf_free( gfld )

      END IF        ! First GRIB2 message

      IF (verbose) WRITE(6,FMT='(1x,I4,2x,I10,1x,I3,3x)',ADVANCE='NO')  &
                   icount,lskip+1,numfields

      DO n = 1, numfields   ! unpack GRIB2 fields in each message

        CALL gf_getfld(cgrib,lgrib,n,.TRUE.,.TRUE.,gfld,iret)
        IF (iret /= 0) THEN
          WRITE(6,'(a,I3)') 'ERROR: Decoding GRIB2 message = ',iret
          GOTO 9999
!          CALL arpsstop('ERROR in rdnmcgrb2 when calling gf_getfld.',1)
        END IF

        IF (verbose) THEN
          IF (n == 1) THEN
          WRITE(6,FMT='(I5,6x,3(a,I3),a,I10,4x)',ADVANCE='NO') n,         &
          '(',gfld%discipline,',',gfld%ipdtmpl(1),',',gfld%ipdtmpl(2),')',&
          gfld%ipdtmpl(10)
          ELSE
          WRITE(6,FMT='(24x,I5,6x,3(a,I3),a,I10,4x)',ADVANCE='NO') n,     &
          '(',gfld%discipline,',',gfld%ipdtmpl(1),',',gfld%ipdtmpl(2),')',&
          gfld%ipdtmpl(10)
          END IF
        END IF

        DO m = 1,maxvar    ! match the required variables
          IF (gfld%discipline == varids(1,m) .AND.      &  ! Discipline
              gfld%ipdtmpl(1) == varids(2,m) .AND.      &  ! Category
              gfld%ipdtmpl(2) == varids(3,m) .AND.      &  ! Parameter
              gfld%ipdtmpl(10)== varids(4,m) ) THEN        ! Elevation

            levelin = gfld%ipdtmpl(12)

!-----------------------------------------------------------------------
!
! Check 3D levels
!
!-----------------------------------------------------------------------

            IF (var3dindx(m) > 0 .AND.         &  ! 3D level
                (gfld%ipdtmpl(10) == 100 .OR.  &  ! Pressure level
                 gfld%ipdtmpl(10) == 105 .OR.  &  ! Hybrid level
                 gfld%ipdtmpl(10) == 104)      &  ! Sigma level
                ) THEN
              DO k = 1,nz_ext
                IF (levelin == var3dlvl(k)) EXIT
              END DO
              IF (k > nz_ext) THEN
                WRITE(6,'(/,1x,a,F15.2,a,/,1x,a,/)')                      &
                'WARNING: variable level (',levelin,') is not found in var3dlvl.', &
                '         Please check the parameters array passing in to rdnmcgrb2.'
                EXIT
  !              CALL arpsstop('Unknown variable level.',1)
              ELSE
                IF (verbose) WRITE(6,FMT='(F10.2,5x,a,I3)',ADVANCE='NO')  &
                                  levelin, '3D-',k
              END IF
              !ij = 0
              DO j = jbgn, jend, jinc
                DO i = ibgn, iend, iinc
                  ij = ( j*jinc+joff )*nx_grb+( i*iinc+ioff )
                  ii = i-imin
                  jj = j-jmin
                  var_grb3d(ii,jj,k,var3dindx(m)) = gfld%fld(ij)  ! store the unpacked 2D slab
                END DO
              END DO
!-----------------------------------------------------------------------
!
! Check Soil levels
!
!-----------------------------------------------------------------------

            ELSE IF (gfld%ipdtmpl(10) == 106 ) THEN  ! Depth below land surface

              DO k = 1,nzsoilin_ext
                IF (levelin == var3dsoil(k)) EXIT
              END DO
              IF (k > nzsoilin_ext) THEN
                IF (lvldbg > 90) WRITE(6,'(/,1x,a,F15.2,a,/,1x,a,/)')                &
                'WARNING: variable level (',levelin,') is not found in var3dsoil.',  &
                '         Please check the parameters array passing in to rdnmcgrb2.'
                EXIT
  !              CALL arpsstop('Unknown soil variable level.',1)
              ELSE
                IF (verbose) WRITE(6,FMT='(F10.2,5x,a,I3.3)',ADVANCE='NO')  &
                                  levelin, '3DSOIL-',k
              END IF
              !ij = 0                            ! assume ijscan = 0
              DO j = jbgn, jend, jinc
                DO i = ibgn, iend, iinc
                  !ij = ij+1
                  ij = ( j*jinc+joff )*nx_grb+( i*iinc+ioff )
                  ii = i-imin
                  jj = j-jmin
                  var_grb3d(ii,jj,k,var3dindx(m)) = gfld%fld(ij)  ! store the unpacked 2D slab
                END DO
              END DO

!-----------------------------------------------------------------------
!
! Check 2D level, Should be at specified 2D slab
!
!-----------------------------------------------------------------------

            ELSE IF (var2dindx(m) > 0 .AND.             & ! 2D level
                     (gfld%ipdtmpl(10) == 103 .OR.      &  ! Specified height above ground
                      gfld%ipdtmpl(10) == 104 )         &  ! Sigma Level
              ) THEN

              IF (levelin == var2dlvl(var2dindx(m)) ) THEN
                !ij = 0
                DO j = jbgn, jend, jinc
                  DO i = ibgn, iend, iinc
                    !ij = ij+1
                    ij = ( j*jinc+joff )*nx_grb+( i*iinc+ioff )
                    ii = i-imin
                    jj = j-jmin
                    var_grb2d(ii,jj,var2dindx(m)) = gfld%fld(ij)
                  END DO
                END DO
                IF (verbose) WRITE(6,FMT='(F10.2,5x,a)',ADVANCE='NO')     &
                                  levelin, '2D'
              ELSE
                IF (verbose) WRITE(6,FMT='(F10.2,5x,a)',ADVANCE='NO')     &
                                  levelin, '-- Skipped'
              END IF  ! 2D variable at the right layer

!-----------------------------------------------------------------------
!
! Check 2D level, Do not need to check the level values
!
!-----------------------------------------------------------------------

            ELSE IF (gfld%ipdtmpl(10) == 1 .OR.             &  ! Ground or water surface
                     gfld%ipdtmpl(10) == 101 ) THEN            ! Mean Sea Level
              !ij = 0
              DO j = jbgn, jend, jinc
                DO i = ibgn, iend, iinc
                  !ij = ij+1
                  ij = ( j*jinc+joff )*nx_grb+( i*iinc+ioff )
                  ii = i-imin
                  jj = j-jmin
                  var_grb2d(ii,jj,var2dindx(m)) = gfld%fld(ij)
                END DO
              END DO
              IF (verbose) WRITE(6,FMT='(F10.2,5x,a,I2,a)',ADVANCE='NO')  &
                                levelin, '2D-ground(',var2dindx(m),')'
            ELSE
              IF (verbose) WRITE(6,FMT='(F10.2,5x,a)',ADVANCE='NO')       &
                                levelin, '-- Skipped'
            END IF         ! Store if matched levels also

            EXIT           ! already matched, so exit from the required variable loop

          END IF         ! matched var. ids
        END DO         ! matching loop all required variables

        IF (verbose ) THEN
          IF ( m > maxvar) THEN
            WRITE(6,'(10x,a)') ' --- skipped.'
          ELSE
            WRITE(6,*)
          END IF
        END IF

        CALL gf_free( gfld )

      END DO               ! numfields

    END DO   ! number of GRIB2 messages

!-----------------------------------------------------------------------
!
! Close file and release the gribfield before returning
!
!-----------------------------------------------------------------------

    CALL baclose( grbunit,iret )
    CALL retunit( grbunit )

!    DEALLOCATE(cgrib)

    IF ( .NOT. ALLOCATED(cgrib) ) THEN
      WRITE(6,'(1x,a,/)') 'WARNING:  File is empty!'
      iret = -888
    ELSE
      DEALLOCATE(cgrib)
    END IF

  END IF                  ! IF (myproc == 0)

  9999  CONTINUE

  CALL mpupdatei(iproj_grb,1)
  CALL mpupdatei(nx_grb,1)
  CALL mpupdatei(ny_grb,1)
  CALL mpupdatei(dx_grb,1)
  CALL mpupdatei(dy_grb,1)
  CALL mpupdater(latsw,1)
  CALL mpupdater(lonsw,1)
  CALL mpupdater(lattru1,1)
  CALL mpupdater(lattru2,1)
  CALL mpupdater(lontrue,1)
  CALL mpupdatei(uvearth,1)
  CALL mpupdater(var_grb2d,nx_ext*ny_ext*n2dvs)
  CALL mpupdater(var_grb3d,nx_ext*ny_ext*nz_ext*n3dvs)
  CALL mpupdatei(iret,1)

  RETURN
END SUBROUTINE rdnmcgrb2
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDGRB2DIMS                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdgrb2dims(gribfile,grbflen,nx_grb,ny_grb,iret)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine extracts dimension sizes from GRIB2 file
!
!  The decoder of GRIB2 is from NCEP.
!
!  NOTE: It is intended to be called from root processor only.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/20/2007
!
!  MODIFICATIONS:
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
!    nx_grb        x-direction grid length
!    ny_grb        y-direction grid length
!
!    iret          Return flag
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  USE GRIB_MOD

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: grbflen

  CHARACTER(LEN=*), INTENT(IN) :: gribfile

  INTEGER, INTENT(OUT) :: nx_grb       ! Number of points along x-axis
  INTEGER, INTENT(OUT) :: ny_grb       ! Number of points along y-axis

  INTEGER, INTENT(OUT) :: iret         ! Return flag

!
!-----------------------------------------------------------------------
!
!  Temporary arrays to read GRIB file
!
!-----------------------------------------------------------------------
!
  INTEGER :: grbunit

  INTEGER, PARAMETER :: maxlen = 32000
  CHARACTER(LEN=1), ALLOCATABLE :: cgrib(:)

  TYPE(gribfield) :: gfld

!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
!  LOGICAL             :: verbose = .FALSE.
  CHARACTER(LEN=12), PARAMETER :: cmd     = 'rdgrib2dims: '

  INTEGER :: iseek
  INTEGER :: lgrib, lgribin, lskip

  INTEGER :: iscan, jscan, ijscan
  INTEGER :: listsec0(3), listsec1(13)
  INTEGER :: numfields, numlocal, maxlocal

  INTEGER :: itemp
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Open the GRIB file
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(1x,8A)') cmd, 'Opening GRIB ', gribfile(1:grbflen)

  CALL getunit( grbunit )

  CALL baopenr(grbunit,gribfile(1:grbflen),iret)

  IF  ( iret /= 0 )  THEN
    WRITE(6,'(a,a,/,a,i4)')                                             &
        'ERROR in opening file ',gribfile(1:grbflen),                   &
        'function baopenr return status: ',iret
    RETURN
  END IF

!-----------------------------------------------------------------------
!
! Loop over all records
!
!-----------------------------------------------------------------------

  iseek   = 0
  iret = 0

  CALL skgb(grbunit,iseek,maxlen,lskip,lgrib)
  IF (lgrib == 0) THEN
    iret = -3
    RETURN
  END IF

  ALLOCATE(cgrib(lgrib), STAT = iret)

  CALL baread(grbunit,lskip,lgrib,lgribin,cgrib)   ! read a GRIB message
  IF (lgrib /= lgribin) THEN
    WRITE(6,'(2(a,I6))')                                              &
       'ERROR: GRIB2 message size in error, required size = ',lgrib,  &
       ', read-in size = ',lgribin
    iret = -3
    RETURN
!    CALL arpsstop('ERROR in rdgrb2dims when calling baread.',1)
  END IF

  CALL gb_info(cgrib,lgrib,listsec0,listsec1,                         &
               numfields,numlocal,maxlocal,iret)
  IF (iret /= 0) THEN
    WRITE(6,'(a,I3)') 'ERROR: quering GRIB2 message = ',iret
    RETURN
!    CALL arpsstop('ERROR in rdgrb2dims when calling gb_info.',1)
  END IF

  CALL gf_getfld(cgrib,lgrib,1,.TRUE.,.TRUE.,gfld,iret)
  IF (iret /= 0) THEN
    WRITE(6,'(a,I3)') 'ERROR: Decoding GRIB2 message = ',iret
    RETURN
!    CALL arpsstop('ERROR in rdgrb2dims when calling gf_getfld.',1)
  END IF

  itemp = gfld%ipdtmpl(5)         ! indicator of model
  IF (gfld%idsect(1) == 7) THEN
    IF (itemp == 83 .OR. itemp == 84) THEN
      WRITE(6,'(1x,a,/)') 'NCEP NAM model'
    ELSE IF (itemp == 81 .OR. itemp == 96) THEN
      WRITE(6,'(1x,a,/)') 'NCEP GFS model'
    ELSE
      WRITE(6,'(1x,a,I3,/)') 'Unkown model from NCEP, model = ',itemp
    END IF
  ELSE
    WRITE(6,'(1x,a,I2,a,I2,/)')                                   &
          'Unknown model and orig center, center = ',gfld%idsect(1),  &
          ' model = ',itemp
  END IF

  IF (gfld%igdtnum == 0) THEN       ! Lat/Lon grid
!    iproj_grb = 4
    nx_grb    = gfld%igdtmpl(8)
    ny_grb    = gfld%igdtmpl(9)
!    dx_grb    = gfld%igdtmpl(17)/1.0E6  ! in degree
!    dy_grb    = gfld%igdtmpl(18)/1.0E6  ! in degree
!    iscan     = IAND(gfld%igdtmpl(19),128)/128  ! get the first bit
                    ! = 0, +x(+i) direction, = 1, -x(-i) direction
!    jscan     = IAND(gfld%igdtmpl(19),64) / 64   ! get the second bit
                    ! = 0, -y(-j) direction, = 1, +y(+j) direction
!    ijscan    = IAND(gfld%igdtmpl(19),32) / 32   ! get the third bit
                    ! = 0, x direction first, = 1, y direction first

!    uvearth   = IAND(gfld%igdtmpl(14),8)  /  8   ! get the fifth bit, big-endian
!    IF (jscan == 1) THEN          ! +y(+j) direction
!      latsw     = gfld%igdtmpl(12)/1.0E6
!      lonsw     = gfld%igdtmpl(13)/1.0E6
!    ELSE                          ! -y(-j) direction
!      latsw     = gfld%igdtmpl(15)/1.0E6
!      lonsw     = gfld%igdtmpl(13)/1.0E6
!    END IF
  ELSE IF (gfld%igdtnum == 30) THEN ! Lambert Conformal Grid
!    iproj_grb   = 2
    nx_grb      = gfld%igdtmpl(8)
    ny_grb      = gfld%igdtmpl(9)
!    lontrue     = gfld%igdtmpl(14)
!    lattru1     = gfld%igdtmpl(19)
!    lattru2     = gfld%igdtmpl(20)
!    dx_grb      = gfld%igdtmpl(15)
!    dy_grb      = gfld%igdtmpl(16)
!    latsw       = gfld%igdtmpl(10)
!    lonsw       = gfld%igdtmpl(11)

  ELSE IF (gfld%igdtnum == 20) THEN ! Polar-Stereographic Grid
!    iproj_grb   = 1
    nx_grb      = gfld%igdtmpl(8)
    ny_grb      = gfld%igdtmpl(9)
!    lattru1     = 60.
!    lattru2     = 91.
!    latsw       = gfld%igdtmpl(10)
!    lonsw       = gfld%igdtmpl(11)

  ELSE
    WRITE(6,'(1x,a,I3)') 'Unkown projection: ',gfld%igdtnum
    iret = -3
    RETURN
!     CALL arpsstop('Unknown map project in GRIB2 file.',1)
  END IF

  CALL gf_free( gfld )

!-----------------------------------------------------------------------
!
! Close file and release the gribfield before returning
!
!-----------------------------------------------------------------------

  CALL baclose( grbunit,iret )
  CALL retunit( grbunit )

  DEALLOCATE (cgrib)

  RETURN
END SUBROUTINE rdgrb2dims
