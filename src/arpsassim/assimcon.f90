!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASSIMCON                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE assimcon(nx,ny,nz,x,y,z,zp,                                  &
           ubar,vbar,pbar,ptbar,rhostr,qvbar,                           &
           u,v,w,ptprt,pprt,qv,qc,qr,qi,qs,qh,                          &
           j1,j2,j3,                                                    &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                     &
           tem9,tem10,tem11)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinates the input of control parameter options for a thermo-
!  dynamic recovery and variational adjustment. The following is a
!  summary of the available options and a list of possible combinations.
!
!  recovopt   Dynamic retrieval option:
!
!          = 0, Do NOT retrieve p' or T'
!          = 1, Retrieve p' and T'       (Implied insertion of the
!                                         interpolated velocities at the
!                                         recovery time.)
!
!  varopt  Variational adjustment option:
!
!          = 0, NO variational adjustment
!          = 1, Perform variational adjustment
!
!  insrtopt   Insertion option:
!
!          = 0, Do NOT insert velocities (NOTE: velocities will be inserted
!                                         in conjunction with dynamic
!                                         retrieval if recovopt=1)
!          = 1, Insert velocities        (At all available data times -
!                                         See NOTE below.)
!
!  A matrix describing these assimilation options is as follows:
!
!    recovopt varopt insrtopt  Comments
!    ----------------------------------------------------------------------
!        0       0       0        Nothing - Exit assimilation mode
!        0       0       1        Direct insertion Only
!        0       1       0        NOT ALLOWED - Exit assimilation mode
!        0       1       1        Var. adj. + direct insertion (no recovery)
!        1       0       1        Recovery + direct insertion (no var. adj.)
!        1       1       1        Recovery, adj., & direct insertion
!        1       1       0        Recovery + var. adj. (no direct insertion)
!        1       0       0        Recovery only
!
!  NOTE:  The insertion option (insrtopt=1) implies that the velocities
!        are to be inserted directly into the model stream at the model
!       time corresponding to that of the data. This inserted velocity
!       data may (varopt=1) or may not (varopt=0) be 'adjusted'.
!       When a recovery is performed (recovopt=1), the 3 input data
!       files (at times t1,t2,& t3) are interpolated to model time
!       steps m1,m2,& m3. In effect, this is also an insertion since
!       the model's prognostic variables at these times are overwritten.
!       Thus, even if the insertion option is not chosen (i.e., insrtopt=0)
!       there will still be an 'insertion' at the recovery time if
!       the recovery option is turned on (recovopt=1).
!
!                  t1             t2             t3
!                 --x-----------|--x--|-----------x--
!                              m1 m2 m3
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Steven Lazarus and Alan Shapiro
!
!  2/24/1993.
!
!  MODIFICATION HISTORY:
!
!  03/08/96 (Limin Zhao)
!  Added an adjustment to the location of radar station so that it
!  will never be on grid points.
!
!  03/19/96 (Limin Zhao)
!  Fixed a bug on assimtim. They need be saved after they are read
!  in first time.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!  OUTPUT:
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    pbar     Base state pressure (Pascal)
!    ptbar    Base state potential temperature (K)
!    rhostr   Base state density (kg/m**3) * j3
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    u        x velocity component at 3 time levels
!    v        y velocity component at 3 time levels
!    w        z velocity component at 3 time levels
!
!    ptprt    Perturbation potential temperature at 3 time levels
!    pprt     Perturbation pressure at 3 time levels
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rain water mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!  NOTE:  Since ASSIMCON reads input data file headers only, the
!       following are dummy arrays.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nt                ! The no. of t-levels of t-dependent arrays

  PARAMETER (nt=3)

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: ubar  (nx,ny,nz)     ! Base state x-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state y-velocity (m/s)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3) times j3
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity (kg/kg)

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qc    (nx,ny,nz,nt)  ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz,nt)  ! Rain water mixing ratio (kg/kg)
  REAL :: qi    (nx,ny,nz,nt)  ! Cloud ice mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz,nt)  ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz,nt)  ! Hail mixing ratio (kg/kg)

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x )
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y )
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z )

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem10  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem11  (nx,ny,nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!

  INTEGER :: ichk     ! Flag indicating whether the next data set is
  SAVE ichk        ! available (0=not available,1=available)
  DATA ichk /1/

  INTEGER :: iassim   ! Flag indicating whether model time (curtim)
  SAVE iassim      ! exceeds the input data time (assimtim).
  DATA iassim /1/  ! (0=yes, exit assimilation mode)

  INTEGER :: icount   ! Flag controling the reading of the batch file.
  SAVE icount      ! (0=do not read batch file,1=read batch file).
  DATA icount /1/

  INTEGER :: iread    ! Counter indicating the input data file and data
  SAVE iread       ! time. Used to flag whether a variational adj.
  DATA iread /1/   ! is to be performed this time step.

  INTEGER :: igo      ! Counter indicating the input data file and data
  SAVE igo         ! time. Used to flag whether a direct insertion
  DATA igo /1/     ! is to be performed this time step.

  INTEGER :: dtcount
  SAVE dtcount
  DATA dtcount /0/

  INTEGER :: isrc     ! Flag coordinating the data ingest
  INTEGER :: istat    ! Status of input data file

  INTEGER :: itag     ! If test flag (itag=ii,igo, or iread)

  REAL :: assimtim(100)     ! Time of (all) input data files
  SAVE assimtim

  REAL :: t1           ! Time of first data file
  REAL :: t2           ! Time of second data file
  REAL :: t3           ! Time of third data file

  REAL :: temscl       ! Grid scale used to calculate cdvdmp

  REAL :: radloc

  REAL :: stor1,stor2,stor3,stor4
  INTEGER :: stor5
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model execution
  INCLUDE 'assim.inc'       ! Recovery/assimilation parameters
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
!  external assimbat
  EXTERNAL assimrd
  EXTERNAL retrint
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
!  If iassim is zero, then one of four things have occured:
!    1) The next available data file to be used in the assimilation
!       arrived late (i.e., after the current model time exceeds that
!       of the input data).  Exit assimilation mode.
!    2) The wrong combination of options for recovopt,velocity,and
!       insrtopt were chosen (see below). Exit assimilation mode.
!    3) The number of files (nvf) has been exceeded.
!    4) The assimilation mode has been turned off (recovopt,varopt,
!       instopt= 0).
!
!-----------------------------------------------------------------------
!
  WRITE(6,*)'code in assimcon, and curtim is   ',curtim

  IF ( iassim == 0) RETURN
!
!-----------------------------------------------------------------------
!
!  Enter this IF block if this is the first entry into ASSIMCON
!  (icount=1).
!
!-----------------------------------------------------------------------
!

  IF ( icount == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Read in user-supplied recovery/assimilation parameters (ASSIMBAT.F)
!  upon first entry into ASSIMCON only (icount=1). The parameter
!  icount is then set to zero to avoid reading the batch file on
!  subsequent visits to this routine.
!
!-----------------------------------------------------------------------
!
!    CALL initassim(nx,ny,nz)

    icount = 0

!
!-----------------------------------------------------------------------
!
!  Check if the radar station located on grid point. If so, shift
!  a small distance. Here we choose to shift (grid space/20).
!
!-----------------------------------------------------------------------
!
!   radloc = xshift - dx*int(xshift/dx)
!   IF (radloc.eq.0) THEN

    radloc = xshift - dx*nint(xshift/dx)
    IF (radloc < dx*1.0E-3) THEN
      xshift = xshift + 0.05*dx
      WRITE(6,*)'radar locates at x grid point'
    END IF

!   radloc = yshift - dy*int(yshift/dy)
!   IF (radloc.eq.0) THEN

    radloc = yshift - dy*nint(yshift/dy)
    IF (radloc < dy*1.0E-3) THEN
      yshift = yshift + 0.05*dy
      WRITE(6,*)'radar locates at y grid point'
    END IF

!   radloc = zshift - dz*int(zshift/dz)
!   IF (radloc.eq.0) THEN

    radloc = zshift - dz*nint(zshift/dz)
    IF (radloc < dz*1.0E-3) THEN
      zshift = zshift + 0.05*dz
      WRITE(6,*)'radar locates at z grid point'
    END IF

    WRITE(6,*)'adjusted radar locations xshift,yshift,zshift=',         &
                                        xshift,yshift,zshift
!
!-----------------------------------------------------------------------
!
!  Check flags to see what combination of assimilation options
!  have been chosen. If an unavailable option is chosen or the
!  assimilation mode has been turned off print message and exit.
!  (See table in Purpose Section above)
!
!-----------------------------------------------------------------------
!
    IF ( recovopt == 0.AND.insrtopt == 0) GO TO 100
!
!-----------------------------------------------------------------------
!
!  Reset nstep so that the I/O corresponds to the restart time and not
!  t=0. This is a temporary fix until the option is introduced into
!  ARPS proper.
!
!-----------------------------------------------------------------------
!
    nstep  = INT( ((curtim-tstart)+0.1*dtbig)/dtbig ) + 1
    PRINT *,'in ASSIMCON, nstep =',nstep
!
!-----------------------------------------------------------------------
!
!  Initialize the counter, ii, to 3. ii is declared in the common
!  block of ASSIM.INC and is used by this routine and ASSIMRD.F only.
!  NOTE:  ii is updated in ASSIMCON only.
!
!-----------------------------------------------------------------------
!
    ii = 3
!
!-----------------------------------------------------------------------
!
!  Read three input data file headers. This is done so that the
!  recovery can access the times of these data files.  The time
!  for each file is brought back via the array 'assimtim' which
!  is then used to set flags to determine when to perform a
!  variational adjustment, and/or direct insertion and/or the first
!  recovery.
!
!  NOTE:  This package assumes that AT LEAST 3 input data files exist
!         prior to model execution.
!
!         The flag isrc coordinates the ingestion of data, where:
!
!         isrc = 1, Flags ASSIMRD to read 3 file headers only
!         isrc = 2, Flags ASSIMRD to read 1 file header  only
!         isrc = 3, Flags ASSIMRD to read 1 complete data file
!                   (called in ASSIMVEL.F)
!         isrc = 4, Flags ASSIMRD to read 3 complete data files
!
!  Note: tem11 is used to read rhobar
!
!-----------------------------------------------------------------------
!
    isrc = 1  ! Flags assimrd to read 3 file headers only
    itag = ii ! Counter indicating the input file name and time

    PRINT *, 'ASSIMCON: CALL assimrd_1 ', itag
    CALL assimrd(nx,ny,nz,x,y,z,zp,                                     &
                 isrc,itag,istat,assimtim,                              &
                 ubar,vbar,pbar,ptbar,rhostr,qvbar,tem11,               &
                 u,v,w,pprt,ptprt,qv,qc,qr,qi,qs,qh,j1,j2,j3,           &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,               &
                 tem9,tem10)

    IF (istat /= 0) THEN
      iassim = 0
      WRITE(6,'(/5x,a,/5x,a,/5x,a,/5x,a)')                              &
          'WARNING: The assimilation package assumes that 3 input',     &
          'data files exist prior to the start of the assimilation.',   &
          'Please provide at least 3 files and restart assimilation.',  &
          'Exit assimilation mode.'
      RETURN
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  If ichk=0, then the data file needed for the next 'assimilation'
!  does not yet exist and thus we set irecov, ivar, and insrt to zero
!  during that model time step. This is done since irecov, ivar, &
!  insrt cannot be flagged using an 'unknown' time (assimtim).
!  If the data file 'arrives' late (i.e., after the model time
!  exceeds that of the input data file), the parameter iassim is set
!  to 0 and the assimilation mode is turned off for the remainder of
!  the simulation. If the data arrive on time, ichk is set to 1, and
!  we proceed below to determine whether a recovery/assimilation is
!  to be performed during this model time step.
!
!-----------------------------------------------------------------------
!
!
!   assimtim(1) =   0.0           ! for 18:36, 18:41, 18:47Z
!   assimtim(2) = 350.0
!   assimtim(3) = 701.0

!  assimtim(1) =  30.0            ! for SSW NOV 1998 tests
!  assimtim(2) = 378.0          ! comment out 12-5-98 SSW
!  assimtim(3) = 732.0

  WRITE(6,*)'assim time are:',assimtim(1),assimtim(2),assimtim(3)

  IF ( ichk == 0) THEN

    irecov  = 0
    ivar = 0
    insrt  = 0

    IF ( ii > nvf.OR.igo > nvf.OR.iread > nvf) iassim = 0

  ELSE

!
!-----------------------------------------------------------------------
!
!  If recovopt equals 1, set flags to perform a retrieval.
!
!-----------------------------------------------------------------------
!
    IF ( recovopt == 1) THEN
!
!-----------------------------------------------------------------------
!
!    Determine whether a recovery/assimilation is to be performed
!    during this model time step.
!
!-----------------------------------------------------------------------
!
      WRITE(6,*) 'curtim = ',curtim
      WRITE(6,*) 'assim time = ',assimtim(1),assimtim(2),assimtim(3)


      IF (curtim < assimtim(ii).OR.ii-1 >= nvf) THEN

        irecov = 0 ! Flag tells RETRPTPR.F not to perform a recovery

      ELSE

        ii = ii + 1
        irecov = 1 ! Flag tells RETRPTPR.F to perform a recovery

      END IF

    END IF
!
!-----------------------------------------------------------------------
!
!  If varopt equals 1, set flags to perform a variational adjustment.
!
!-----------------------------------------------------------------------
!

    IF (varopt == 1) THEN

!
!-----------------------------------------------------------------------
!
!    Determine whether a variational adjustment is to be performed
!    during this model time step.
!
!-----------------------------------------------------------------------
!
      IF (curtim < assimtim(iread) .OR. iread-1 >= nvf) THEN

        ivar = 0 ! Flag tells ASSIMVEL.F not to perform a variational adjustment

      ELSE

        iread = iread + 1

        ivar = 1 ! Flag tells ASSIMVEL.F to perform a variational adjustment

      END IF

    END IF

!
!-----------------------------------------------------------------------
!
!  If insrtopt equals 1, set flags to perform a direct insertion.
!
!-----------------------------------------------------------------------
!

    IF (insrtopt == 1) THEN

!
!-----------------------------------------------------------------------
!
!    Determine whether a direct insertion is to be performed
!    during this model time step.
!
!-----------------------------------------------------------------------
!
      IF (curtim < assimtim(igo) .OR. igo-1 >= nvf) THEN

        insrt = 0 ! Flag tells ASSIMVEL.F not to perform an insertion

      ELSE

        igo = igo + 1

        insrt = 1 ! Flag tells ASSIMVEL.F to perform an insertion

      END IF

    END IF

  END IF

  WRITE(6,*) 'irecov,ivar,insrt: ',irecov,ivar,insrt
  WRITE(6,*) 'ii,iread,igo: ',ii,iread,igo
!
!-----------------------------------------------------------------------
!
!    Read the next available data file. This is done in order to set the
!    flags, irecov, ivar, and insrt for the next recovery, variational
!    adjustment, and/or insertion.  Since it is assumed that at least 3
!    data files exist prior to the assimilation period, we do not read
!    the 4th file header until the model time exceeds the time of the
!    3rd file (t3).
!
!          t1             t2             t3             t4
!         --x--------------x--------------x--------------x--
!                                          |
!                               Begin search for 4th File
!
!    Similarly, we do not begin looking for the next input data file
!    header [at time, t(n)] until the model time exceeds the time of the
!    preceeding input data file [at time, t(n-1)].
!
!                        t(n-1)         t(n)
!                        --x--------------x--
!                           |
!                Begin search for nth File
!
!-----------------------------------------------------------------------
!
  IF ( ii > 3.OR.iread > 3.OR.igo > 3) THEN

    IF ( irecov == 1.OR.ivar == 1.OR.insrt == 1.OR.ichk == 0) THEN

      isrc   = 2 ! Flags assimrd to read 1 file header only

      itag   = MAX(ii,iread,igo)

      PRINT *, 'ASSIMCON: CALL assimrd_2 ',                             &
          itag, irecov, ivar, insrt
      CALL assimrd(nx,ny,nz,x,y,z,zp,                                   &
                   isrc,itag,istat,assimtim,                            &
                   ubar,vbar,pbar,ptbar,rhostr,qvbar,tem11,             &
                   u,v,w,pprt,ptprt,qv,qc,qr,qi,qs,qh,j1,j2,j3,         &
                   tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,             &
                   tem9,tem10)


      IF (istat /= 0 .OR. itag > nvf) THEN

        ichk   = 0

      ELSE

        ichk = 1
!
!-----------------------------------------------------------------------
!
!  Check to ensure that the input of the next data file occurs no more
!  than one model time step (dtbig) past the next recovery time
!  assimtim(ii). If it does, exit assimilation mode by setting
!  iassim = 0.
!  NOTE:  This precludes any further assimilation even if subsequent
!       data files become available.
!
!-----------------------------------------------------------------------
!
        IF (curtim-assimtim(itag) > dtbig) THEN

          iassim = 0

          irecov  = 0
          ivar = 0
          insrt  = 0

          WRITE(6,'(/5x,a,/5x,a,/5x,a)')                                &
              'WARNING: The successful read of the input data file occured', &
              'after the insertion/adjustment/recovery window',         &
              'Exit assimilation mode.'

          RETURN

        END IF

      END IF

    END IF

  END IF

  RETURN

  100   WRITE(6,'(/5x,a,/5x,a,/10x,a,i2,/10x,a,i2,/10x,a,i2,/5x,a,/5x,a)') &
        'WARNING: You have either chosen an option that is not available', &
        'or have turned the assimilation mode off. The input option was', &
        'recovopt=',recovopt,'varopt=',varopt,'insrtopt=',insrtopt,     &
        'See ASSIMCON for appropriate parameters.',                     &
        'Exit assimilation mode.'

  iassim = 0

  RETURN

END SUBROUTINE assimcon
