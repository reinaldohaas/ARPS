
PROGRAM split_general
!-----------------------------------------------------------------------
!
! PURPOSE: 
!    Split any files in either HDF 4 format or netCDF format into
!    patches. The files can contain 1d, 2d, 3d or 4d either integer
!    arrays or float arrays. The patched files will contain the same
!    data as original file but in evenly divided subdomain specified
!    by the user.
!
!-----------------------------------------------------------------------
!
! Author: Yunheng Wang (10/27/2006)
!
! MODIFICATIONS:
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE
!-----------------------------------------------------------------------
!
! NAMELIST variables
!
!-----------------------------------------------------------------------

  INTEGER :: finfmt   ! = 3, for HDF4
                      ! = 7, for netCDF
  INTEGER :: finopt   ! = 1, regular files
                      ! = 2, explicit list

  CHARACTER(LEN=256) :: fheader         ! for finopt = 1 only
  CHARACTER(LEN=256) :: ftrailer
  REAL    :: tintv_in
  REAL    :: tbgn_in
  REAL    :: tend_in

  INTEGER, PARAMETER :: nfile_max = 100 ! for finopt = 2
  INTEGER            :: nfile           ! also used for general purpose 
  CHARACTER(LEN=256) :: filenames(nfile_max)

  NAMELIST /file_names/ finfmt, finopt, fheader, ftrailer,              &
                        tbgn_in,tintv_in,tend_in, nfile, filenames

  INTEGER :: nproc_x, nproc_y
  LOGICAL :: dimnamein
  CHARACTER(LEN=256) :: xdimname, ydimname
  LOGICAL :: stagdims
  INTEGER :: varidx, nxidx, nyidx
  NAMELIST /message_passing/ nproc_x, nproc_y,                          &
                             dimnamein, xdimname, ydimname, stagdims,   &
                             varidx, nxidx, nyidx

  CHARACTER(LEN=256) :: outdirname
  NAMELIST /output/ outdirname

  INTEGER :: debug
  NAMELIST /debugging/ debug

!
!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nf

  INTEGER :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code below
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  finfmt = 7
  finopt = 2
  fheader  = ' '
  ftrailer = ' '
  tintv_in = 0.0
  tbgn_in  = 0.0
  tend_in  = 0.0
  nfile  = 0
  filenames(:) = ' '

  READ(5,file_names,ERR=999)
  WRITE(6,'(1x,a)'   ) 'Namelist file_names was successfully read.'
  WRITE(6,'(3x,a,i3)') 'Input finopt = ', finopt
  WRITE(6,'(3x,a,i3)') 'Input finfmt = ', finfmt

  IF (finopt == 1) THEN
    CALL getinfns(finfmt, fheader, ftrailer, tintv_in, tbgn_in, tend_in,&
                  filenames, nfile_max, nfile, istatus)
    IF (istatus < 0) STOP
  END IF

  WRITE(6,'(3x,a,i3)') 'Input nfile  = ', nfile

  WRITE(6,'(3x,a)') 'Files to split are:'
  DO nf=1,nfile
    WRITE(6,'(7x,a,i3,a,a)') 'No.',nf,' is ',TRIM(filenames(nf))
  END DO
 
  nproc_x = 1
  nproc_y = 1
  dimnamein = .TRUE.
  xdimname = ' '
  ydimname = ' '
  stagdims = .FALSE.
  varidx = 2
  nxidx  = 1
  nyidx  = 2
  

  READ(5,message_passing,ERR=999)
  WRITE(6,'(/,1x,a)'   ) 'Namelist message_passing was successfully read.'
  WRITE(6,'(3x,a,i3)') 'nproc_x = ', nproc_x
  WRITE(6,'(3x,a,i3)') 'nproc_y = ', nproc_y
  WRITE(6,'(3x,a,l2)') 'dimnamein = ', dimnamein
  WRITE(6,'(3x,a,I2)') '  varidx = ', varidx
  WRITE(6,'(3x,a,I2)') '  nxidx  = ', nxidx
  WRITE(6,'(3x,a,I2)') '  nyidx  = ', nyidx
  WRITE(6,'(3x,a,a)')  'xdimname  = ', TRIM(xdimname)
  WRITE(6,'(3x,a,a)')  'ydimname  = ', TRIM(ydimname)
  WRITE(6,'(3x,a,L2)') 'stagdims  = ', stagdims

  outdirname = './'
  READ(5,output,ERR=999)
  WRITE(6,'(/,1x,a)'   ) 'Namelist output was successfully read.'
  WRITE(6,'(3x,2a)')   'outdirname = ', TRIM(outdirname)

  debug = 0
  READ(5,debugging,ERR=999)
  WRITE(6,'(/,1x,a)'   ) 'Namelist debugging was successfully read.'
  WRITE(6,'(3x,a,i3,/)') 'debug = ', debug

  WRITE(6,'(1x,a,/)') '*****************************'
!-----------------------------------------------------------------------
!
! Calling specific subroutines to do the main job
!
!-----------------------------------------------------------------------

  IF (finfmt == 3) THEN
    CALL splithdf(filenames,nfile,dimnamein,xdimname,ydimname,          &
                  varidx,nxidx,nyidx,nproc_x,nproc_y,outdirname,debug,  &
                  istatus)
  ELSE IF (finfmt == 7 .OR. finfmt == 8) THEN
    CALL splitncdf(filenames,nfile,stagdims,xdimname,ydimname,          &
                   nproc_x,nproc_y,outdirname,debug,                    &
                   istatus)
  ELSE
    WRITE(6,'(1x,a,I2)') 'ERROR: unsupported file format = ',finfmt
    WRITE(6,'(1x,a)')    '       The program only support files '// &
            'in HDF 4 format (finfmt = 3) or netCDF format (finfmt = 7/8)'
    istatus = 1
  END IF

  GO TO 100

!-----------------------------------------------------------------------
!
! Just before termination
!
!-----------------------------------------------------------------------

  999  WRITE(6,'(1x, a,a)') 'Error reading NAMELIST file. Job stopped.'
  STOP

  100    CONTINUE
  IF (istatus == 0) THEN
    WRITE(6,'(/,4x,a,/)') '==== Program SPLIT terminated normally ===='
  ELSE
    WRITE(6,'(/,4x,a,I3,a/)') '**** Program SPLIT terminated with error = ',istatus,' ****'
  END IF

  STOP
END PROGRAM split_general

SUBROUTINE getinfns(finfmt, fheader, ftrailer, tintv_in, tbgn_in, tend_in, &
                    filenames, nfile_max, nfile, istatus)

!-----------------------------------------------------------------------
!
! Purpose:
!   Construct filenames array from input parameters
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: finfmt
  CHARACTER(LEN=*), INTENT(IN)  :: fheader, ftrailer
  REAL,             INTENT(IN)  :: tbgn_in, tintv_in, tend_in
  INTEGER,          INTENT(IN)  :: nfile_max

  CHARACTER(LEN=256), INTENT(OUT) :: filenames(nfile_max)
  INTEGER,            INTENT(OUT) :: nfile
  INTEGER,            INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  CHARACTER(LEN=3) :: fmtstr
 
  INTEGER :: n
  REAL    :: time
  INTEGER :: lheader, ltrailer
  INTEGER :: itime

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (finfmt == 3) THEN
    fmtstr = 'hdf'
  ELSE IF (finfmt == 7 .OR. finfmt == 8) THEN
    fmtstr = 'net'
  ELSE
    WRITE(6,'(/,1x,a,I2,/)') 'ERROR: unsupported data format: ',finfmt
    istatus = -1
    RETURN
  END IF
     
  lheader = LEN_TRIM(fheader)
  ltrailer= LEN_TRIM(ftrailer)

  nfile = 0
  time = tbgn_in
  DO n = 1, nfile_max
    IF (time > tend_in + 0.01*tintv_in) EXIT
    nfile = nfile + 1
    itime = INT(time)
    WRITE(filenames(n),'(3a,I6.6,a)') fheader(1:lheader),'.',fmtstr,    &
                                      itime,ftrailer(1:ltrailer)
    time = tbgn_in + n*tintv_in
  END DO

  RETURN
END SUBROUTINE getinfns
