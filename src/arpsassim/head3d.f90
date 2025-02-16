!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DTAHEAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dtahead(nx,ny,nz,                                            &
           hinfmt ,grdbasfn,lengbf,datafn,lendtf,time,                  &
           x,y,z,zp, uprt ,vprt ,wprt ,ptprt, pprt ,                    &
           qvprt, qc, qr, qi, qs, qh, km,                               &
           ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                &
           ireturn, tem1, tem2, tem3)

!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the reading of history data headers of various formats.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Steven Lazarus
!    2/1/1994.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx,ny,nz The dimension of data arrays
!
!    hinfmt   The format of the history data dump
!             =1, machine dependent unformatted binary dump,
!             =2, formatted ascii dump,
!
!    grdbasfn Name of the grid/base state array file
!    lengbf   Length of the grid/base state data file name string
!    datafn   Name of the other time dependent data file
!    lendtf   Length of the data file name string
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    wprt     vertical component of perturbation velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!
!    qvprt    perturbation water vapor mixing ratio (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    km       Turbulent mixing coefficient (m**2/s)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state air density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!  OUTPUT:
!
!    time     The time of the input data (s)
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    wprt     vertical component of perturbation velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!
!    qvprt    perturbation water vapor mixing ratio (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    km       Turbulent mixing coefficient (m**2/s)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state air density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz
  REAL :: time
  INTEGER :: hinfmt,lengbf,lendtf
  CHARACTER (LEN=1) :: grdbasfn
  CHARACTER (LEN=1) :: datafn
!
  REAL :: x    (nx)
  REAL :: y    (ny)
  REAL :: z    (nz)
  REAL :: zp   (nx,ny,nz)

  REAL :: uprt (nx,ny,nz)
  REAL :: vprt (nx,ny,nz)
  REAL :: wprt (nx,ny,nz)
  REAL :: ptprt(nx,ny,nz)
  REAL :: pprt (nx,ny,nz)
  REAL :: qvprt(nx,ny,nz)
  REAL :: qc   (nx,ny,nz)
  REAL :: qr   (nx,ny,nz)
  REAL :: qi   (nx,ny,nz)
  REAL :: qs   (nx,ny,nz)
  REAL :: qh   (nx,ny,nz)
  REAL :: km    (nx,ny,nz)
  REAL :: ubar  (nx,ny,nz)
  REAL :: vbar  (nx,ny,nz)
  REAL :: wbar  (nx,ny,nz)
  REAL :: ptbar (nx,ny,nz)
  REAL :: rhobar(nx,ny,nz)
  REAL :: pbar  (nx,ny,nz)
  REAL :: qvbar (nx,ny,nz)
!
  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
  REAL :: tem3(nx,ny,nz)

  INTEGER :: ngchan,nchanl,ireturn,istat
  INTEGER :: grdbas

  REAL :: btime ! The time of the base state data
  REAL :: amin, amax
  INTEGER :: i,j,k
  INTEGER :: ierr
  LOGICAL :: fexist
  INTEGER :: packed
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'indtflg.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
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
!  Read data fields.
!
!-----------------------------------------------------------------------
!
  CALL getunit(nchanl)
  grdbas = 0

  INQUIRE(FILE=datafn(1:lendtf), EXIST = fexist )
  IF( fexist ) GO TO 100

  INQUIRE(FILE=datafn(1:lendtf)//'.Z', EXIST = fexist )
  IF( fexist ) THEN
    CALL uncmprs( datafn(1:lendtf)//'.Z' )
    GO TO 100
  END IF

  INQUIRE(FILE=datafn(1:lendtf)//'.gz', EXIST = fexist )
  IF( fexist ) THEN
    CALL uncmprs( datafn(1:lendtf)//'.gz' )
    GO TO 100
  END IF

  WRITE(6,'(/1x,a,/1x,a/)')                                             &
       'File '//datafn(1:lendtf)                                        &
       //' or its compressed version not found.',                       &
       'Program returned from DTAHEAD.'
  RETURN

  100   CONTINUE

  IF( hinfmt == 1 ) THEN

!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(datafn(1:lendtf), '-F f77 -N ieee', ierr)

    OPEN(UNIT=nchanl,FILE=datafn(1:lendtf),                             &
         STATUS='old',FORM='unformatted',IOSTAT=istat)

    IF( istat /= 0 ) GO TO 998

    CALL binhead(nchanl,time,ireturn)

    CLOSE(UNIT=nchanl)

  ELSE IF( hinfmt == 2 ) THEN

    OPEN(UNIT=nchanl,FILE=datafn(1:lendtf),                             &
         STATUS='old',FORM='formatted',IOSTAT=istat)

    IF( istat /= 0 ) GO TO 998

    CALL aschead(nchanl,time,ireturn)

    CLOSE(UNIT=nchanl)

!   ELSEIF( hinfmt .eq.3 ) THEN

!     CALL hdfhead(datafn(1:lendtf),time,ireturn)

!   ELSEIF( hinfmt .eq.4 ) THEN

!     open(unit=nchanl,file=datafn(1:lendtf),
!  :       status='old',form='unformatted',iostat=istat)

!     IF( istat.ne.0 ) GOTO 998

!     CALL pakhead(nchanl,time,ireturn,grdbas)

!     close(unit=nchanl)

!   ELSEIF( hinfmt .eq.6 ) THEN

!     open(unit=nchanl,file=datafn(1:lendtf),
!  :          status='old',form='unformatted',iostat=istat)

!     IF( istat.ne.0 ) GOTO 998

!     CALL bn2head(nchanl,time,ireturn)

!     close(unit=nchanl)

!   ELSE IF (hinfmt .eq. 7) THEN     ! NetCDF format  *NOT AVAILABLE YET*

!     packed = 0
!     CALL nethead (datafn(1:lendtf),time,packed)


!   ELSE IF (hinfmt .eq. 8) THEN     ! NetCDF packed format

!     packed = 1
!     CALL nethead (datafn(1:lendtf),time,packed)

  ELSE

    WRITE(6,'(a,i3,a)')                                                 &
          ' Data format flag had an invalid value ',                    &
           hinfmt ,' program stopped.'
    STOP
  END IF

  CALL retunit(nchanl)

  CLOSE(UNIT=nchanl)

  RETURN

  998   CONTINUE
  WRITE(6,'(1x,a,a,/1x,i3,a)')                                          &
      'Error occured when opening file ',datafn(1:lendtf),              &
      'using FORTRAN unit ',nchanl,' Program returned from DTAHEAD.'

  RETURN
END SUBROUTINE dtahead
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BINHEAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE binhead(inch,time,ireturn)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in data header ONLY from binary dataset created by ARPS.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Steven Lazarus
!  1/28/94.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    inch     Channel number for binary reading.
!             This channel must be opened for unformatted reading
!             by the calling routine.
!
!  OUTPUT:
!
!    time     Time in seconds of data read from "filename"
!
!    ireturn  Return status indicator
!             =0, successful read of all data
!             =1, error reading data
!             =2, end-of-file reached during read attempt
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
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

  REAL :: time

  INTEGER :: inch
  INTEGER :: ireturn
!
!-----------------------------------------------------------------------
!
!  Parameters describing routine that wrote the gridded data
!
!-----------------------------------------------------------------------
!
!
!
!   character*40 fmtver,fmtverin
!   parameter (fmtver='003.20 Binary Data')

  CHARACTER (LEN=40) :: fmtver0,fmtver1,fmtver,fmtverin
  PARAMETER (fmtver='004.10 Binary Data')
  PARAMETER (fmtver0='003.20 Binary Data')
  PARAMETER (fmtver1='004.00 Binary Data')
  INTEGER :: oldver            ! Flag indicating if the file is an old
                               ! (oldver=1) or current format (oldver=0).

  CHARACTER (LEN=10) :: tmunit
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'indtflg.inc'
  INCLUDE 'globcst.inc'
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
!  Read header info
!
!-----------------------------------------------------------------------
!
  READ(inch,ERR=110,END=120) fmtverin

!
!
!   IF( fmtverin .ne. fmtver ) THEN
!     write(6,'(/1x,a,/1x,2a,/1x,3a)')
!  :  'Data format incompatible with the data reader.',
!  :  'Format of data is ',fmtverin,' Format of reader is ',fmtver,
!  :  '. Job stopped.'
!     STOP
!   ENDIF

  IF( (fmtverin /= fmtver) .AND. (fmtverin /= fmtver0)                  &
        .AND. (fmtverin /= fmtver1) ) THEN
    WRITE(6,'(/1x,a/1x,2a/1x,2a/1x,2a/1x,a)')                           &
        'Data format incompatible with the data reader.',               &
        'Format of data is ',fmtverin,' Format of reader is ',fmtver1,  &
        'compitable to ',fmtver0, '. Job stopped.'
    STOP
  END IF

  IF ( fmtverin == fmtver ) THEN
    oldver = 0
  ELSE
    oldver = 1
  END IF

  READ(inch,ERR=110,END=120) runname
  READ(inch,ERR=110,END=120) nocmnt
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      READ(inch,ERR=110,END=120) cmnt(i)
    END DO
  END IF

  WRITE(6,'(//'' THE NAME OF THE INPUT DATA IS:  '',A//)') runname

  READ(inch,ERR=110,END=120) time,tmunit

!
!-----------------------------------------------------------------------
!
!  Exit message
!
!----------------------------------------------------------------------
!
  WRITE(6,'(/a,F8.1,a/)')                                               &
      ' Data header at time=', time,' (sec) were successfully read.'

  ireturn = 0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Error during read
!
!----------------------------------------------------------------------
!

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in BINHEAD'
  ireturn=1
  RETURN
!
!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!----------------------------------------------------------------------
!

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in BINHEAD'
  ireturn=2
  RETURN
END SUBROUTINE binhead
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASCHEAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE aschead(inch,time,ireturn)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read history data header ONLY from channel nchanl in formatted
!  ASCII data format.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Steven Lazarus
!  1/28/94.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    inch     Channel number for ASCII reading.
!             This channel must be opened for formatted reading
!             by the calling routine.
!
!  OUTPUT:
!
!    time     Time in seconds of data read from "filename"
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: inch
  INTEGER :: ireturn

  REAL :: time
!
!-----------------------------------------------------------------------
!
!  Parameters describing routine that wrote the gridded data
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=40) :: fmtver,fmtverin
  PARAMETER (fmtver='003.20 ASCII Formatted Data')
  CHARACTER (LEN=10) :: tmunit
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'indtflg.inc'
  INCLUDE 'globcst.inc'
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
!  Read header info
!
!-----------------------------------------------------------------------
!
  READ(inch,'(1x,a40)',ERR=110,END=120) fmtverin

  IF( fmtverin /= fmtver ) THEN
    WRITE(6,'(/1x,a,/1x,2a,/1x,3a)')                                    &
        'Data format incompatible with the data reader.',               &
        'Format of data is ',fmtverin,' Format of reader is ',fmtver,   &
        '. Job stopped.'
    STOP
  END IF

  READ(inch,'(1x,a80)',ERR=110,END=120) runname
  WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') runname
!
  READ(inch,'(1x,i4)',ERR=110,END=120) nocmnt
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      READ(inch,'(1x,a80)',ERR=110,END=120) cmnt(i)
    END DO
  END IF
!
  READ(inch,'(1x,e16.8,1x,a10)',ERR=110,END=120) time,tmunit

!
!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!----------------------------------------------------------------------
!
  930   CONTINUE

  WRITE(6,'(/a,F8.1,a/)')                                               &
      'Header at time=', time,' (sec) was successfully read.'

  ireturn = 0

  RETURN

!
!-----------------------------------------------------------------------
!
!  Error during read
!
!----------------------------------------------------------------------
!

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in ASCHEAD'
  ireturn=1
  RETURN
!
!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!----------------------------------------------------------------------
!

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in ASCHEAD'
  ireturn=2
  RETURN
END SUBROUTINE aschead
