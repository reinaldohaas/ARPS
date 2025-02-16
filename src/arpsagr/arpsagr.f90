!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPSAGR                    ######
!######             ARPS Adaptive Grid Refinement            ######
!######                                                      ######
!######                     Developed by                     ######
!######           William C. Skamarock (NCAR/MMM)            ######
!######                                                      ######
!######                        AND                           ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

PROGRAM agri
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  3-D ADAPTIVE MODEL INTERFACE
!
!  This is the main driver and includes the top level data read/write
!  routines.  Only the interface routines are in this directory,
!  the solver and user supplied I/O and plotting routines are
!  elsewhere.
!
!  The include files that go along with the interface are
!
!  agrialloc.inc -- common block and size for main storage array.
!  grddsc.inc    -- common block for variable descriptions.
!  agrigrid.inc  -- parameters describing the number of variables,
!                   constants, etc. for any grid.
!  manage.inc    -- storage management common block.
!  nodal.inc     -- common block containing grid location and
!                   relationship information, along with some
!                   adaptive grid timing parameters.
!  agricpu.inc   -- common block containing the timing statistics for
!                   various components of the interface and solver code.
!  agricst.inc   -- Common blocks containing the character strings
!                   for runnames etc.
!                   Common block containing the logicals controlling
!                   diagnostic output.
!                   Common block containing the timing statistics for
!                   various components of the interface and solver code.
!                   Common block to pass parameter defining the
!                   order of spatial interpolations.
!
!  Data definitions for node and rnode are in nodal.inc.
!
!  The input file should be in unit 7.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: William C. Skamarock (NCAR/MMM) and Ming Xue (CAPS/OU)
!
!    For information contact Bill Skamarock at (303) 497-8893
!    or skamaroc@mmm.ucar.edu
!
!
!  MODIFICATIONS:
!
!  10/30/1992 (William C. Skamarock and Ming Xue)
!
!  08/??/1995 (E.J. Adlerman)
!  Used ARPS 4.0.22 as the solver
!
!  03/27/1997 (Yuhe Liu)
!  1. Fully upgraded to ARPS 4.2.4;
!  2. Moved the AGRI input parameters into arps.input by creating
!     a new namelist, &arpsagr and added documentation for the
!     parameters in the new namelist;
!  3. Rewrote the makefile, makefile.agri, in order to be controled
!     by makearps script
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agricpu.inc'
  INCLUDE 'agricst.inc'
!
!-----------------------------------------------------------------------
!
!  Include file for global constants used by ARPS.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'
!
!-----------------------------------------------------------------------
!
!  Define the namelist, &arpsagr, for ARPS AGR
!
!-----------------------------------------------------------------------
!
  NAMELIST /arpsagr/  runold, rstime, nxc, nyc, nzc,                    &
                      levfix, intrat, intratt,                          &
                      intrpodr, kcheck,                                 &
                      verbose1, verbose2, verbose3,                     &
                      verbose4, verbose5, verbose6,                     &
                      rstart,rstdump,grdsrt,                            &
                      nfinelv, ngrdnew,                                 &
                      ixc,jyc,ixln,jyln,gangle

  INTEGER :: iout(6)

  INTEGER :: nstart,nstop
  INTEGER :: i, lv

  INTEGER :: nx,ny,nz       ! Actual dimension sizes
  INTEGER :: ii,ir,mptr
  INTEGER :: igtint,igtrel  ! Functions to return the pointers to the
                            ! integer or real constant array
  REAL :: timeend
  REAL :: tstopnew
  REAL :: f_cputime, dtbase

  LOGICAL :: iexist

!  character*80 cstfile     ! not used
!  logical resetcon         ! not used

  INTEGER :: iorder           ! not used
  REAL :: cut              ! not used
  REAL :: relcut           ! not used
  REAL :: cdist            ! not used
  INTEGER :: nstyps

  CHARACTER(LEN=256) :: namelist_filename
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(/ 16(/5x,a)//)')                                            &
      '###############################################################', &
      '###############################################################', &
      '#####                                                     #####', &
      '#####                      Welcome to                     #####', &
      '#####                                                     #####', &
      '#####      Adaptive Grid Refinement Interface (AGRI)      #####', &
      '#####                                                     #####', &
      '#####                        For                          #####', &
      '#####                                                     #####', &
      '#####   The Advanced Regional Prediction System  (ARPS)   #####', &
      '#####                                                     #####', &
      '#####                     Developed by                    #####', &
      '#####     Center for Analysis and Prediction of Storms    #####', &
      '#####                University of Oklahoma               #####', &
      '#####                                                     #####', &
      '###############################################################', &
      '###############################################################'

  mptr = 1    ! Set the grid number to 1 for base grid.
  mgrid = mptr
  nestgrd = 1 ! Set the grid nesting flag to 1, i.e. nesting.
!
!-----------------------------------------------------------------------
!
!  Call prtcpu(0,6) to initialize the cpu statistics
!
!-----------------------------------------------------------------------
!
  CALL prtcpu(0,6)
  cpu_simulation = f_cputime()
!
!-----------------------------------------------------------------------
!
!  open the graphics stuff for test version
!
!-----------------------------------------------------------------------
!
!  call xdevic
!
!-----------------------------------------------------------------------
!
!  set storage highwater counter
!
!-----------------------------------------------------------------------
!
  ihighwater = 0
!
!-----------------------------------------------------------------------
!
!  set up units for io
!
!-----------------------------------------------------------------------
!
  CALL iosetup
!
!-----------------------------------------------------------------------
!
!  read in headline info and dump back out to standard output
!
!-----------------------------------------------------------------------
!
! READ(5,arpsagr)  ! now read by initpara.

!-----------------------------------------------------------------------
!
!  Read in ARPS model configuration parameters
!
!-----------------------------------------------------------------------
!
  namelist_filename = ' '
  CALL initpara(nxc,nyc,nzc,nstyps, namelist_filename)

  IF( nstyps /= 4 ) THEN
    WRITE (6,'(a,i4,a,i4)') 'ARPS AGR always assumes nstyps <= 4.'
  ENDIF

  nmlntho = 80
  CALL gtlfnkey( runold, nmlntho )

  IF ( nfinelv > nfinelv_max ) THEN
    WRITE (6,'(a,i4,a,i4)')                                             &
        'The number of new grid levels should not be greater ',         &
        'than ',nfinelv_max,', Reset it set to ',nfinelv_max
    nfinelv = nfinelv_max
  ELSE IF ( nfinelv < 0 ) THEN
    WRITE (6,'(a,i4,a,i4)')                                             &
        'The number of new grid levels should not be negative. ',       &
        'Reset it set to ',0
    nfinelv = 0
  END IF

  lfine  = nfinelv + 1
  mxnest = nfinelv + 1

  IF ( nfinelv > 0 ) THEN
    DO lv=1,nfinelv
      IF ( ngrdnew(lv) < 0 ) THEN
        ngrdnew(lv) = 0
        WRITE (6,'(a)')                                                 &
            'The number of new grid should not be less ',               &
            'than 0, Reset it set to 0'
      ELSE IF ( ngrdnew(lv) > ngrdnew_max ) THEN
        ngrdnew(lv) = ngrdnew_max
        WRITE (6,'(a,i4,a,i4)')                                         &
            'The number of new grid should not be greater ',            &
            'than ',ngrdnew_max,', Reset it set to ',ngrdnew_max
      END IF
    END DO
  END IF

  DO lv=nfinelv+1,nfinelv_max
    ngrdnew(lv) = 0
  END DO

! WRITE(6,'(/1x,a)')      '&arpsagr'
! WRITE (6,'(3x,a)')      'runold   = '''//runold(1:nmlntho)//''','
! WRITE (6,'(3x,a,f16.4,a)') 'rstime   = ', rstime,   ','
! WRITE (6,'(3x,a,l3,a)')    'verbose1 = ', verbose1, ','
! WRITE (6,'(3x,a,l3,a)')    'verbose2 = ', verbose2, ','
! WRITE (6,'(3x,a,l3,a)')    'verbose3 = ', verbose3, ','
! WRITE (6,'(3x,a,l3,a)')    'verbose4 = ', verbose4, ','
! WRITE (6,'(3x,a,l3,a)')    'verbose5 = ', verbose5, ','
! WRITE (6,'(3x,a,l3,a)')    'verbose6 = ', verbose6, ','
! WRITE (6,'(3x,a,l3,a)')    'rstart   = ', rstart,   ','
! WRITE (6,'(3x,a,l3,a)')    'rstdump  = ', rstdump,  ','
! WRITE (6,'(3x,a,l3,a)')    'grdsrt   = ', grdsrt,   ','
! WRITE (6,'(3x,a,i4,a)')    'nxc      = ', nxc,      ','
! WRITE (6,'(3x,a,i4,a)')    'nyc      = ', nyc,      ','
! WRITE (6,'(3x,a,i4,a)')    'nzc      = ', nzc,      ','
! WRITE (6,'(3x,a,i4,a)')    'levfix   = ', levfix,   ','
! WRITE (6,'(3x,a,i4,a)')    'intrat   = ', intrat,   ','
! WRITE (6,'(3x,a,i4,a)')    'intratt  = ', intratt,  ','
! WRITE (6,'(3x,a,i8,a)')    'kcheck   = ', kcheck,   ','
! WRITE (6,'(3x,a,i4,a)')    'intrpodr = ', intrpodr, ','
! WRITE (6,'(3x,a,i4,a)')    'nfinelv  = ', nfinelv,  ','

! IF ( nfinelv > 0 ) THEN
!   DO lv=1,nfinelv
!     WRITE (6,'(3x,a,i2.2,a,i4,a)')                                    &
!         'ngrdnew(',lv,') = ', ngrdnew(lv), ','
!     DO i=1,ngrdnew(lv)
!       WRITE (6,'(3x,a,i2.2,a,i2.2,a,f5.1,a,i2.2,a,i2.2,a,f5.1,        &
!       &          3X,a,i2.2,a,i2.2,a,f5.1,a,i2.2,a,i2.2,a,f5.1,        &
!       &          3X,a,i2.2,a,i2.2,a,f5.1,a)')                         &
!           'ixc(',i,',',lv,') = ',ixc(i,lv),                           &
!           ',jyc(',i,',',lv,')=',jyc(i,lv),                            &
!           ',ixln(',i,',',lv,')=',ixln(i,lv),                          &
!           ',jyln(',i,',',lv,')=',jyln(i,lv),                          &
!           ',gangle(',i,',',lv,')= ',gangle(i,lv),','
!     END DO
!   END DO
! END IF
!

  runnew = runname
  nmlnthn = lfnkey
  tstopnew = tstop      ! to prevent overwrite the value from restart

  IF ( hdmpfmt == 5 .OR. hdmpfmt == 9 ) THEN
    WRITE(6,'(a,i2/a)')                                                 &
        ' ARPSAGRI does not support hidtory dump format ',              &
        hdmpfmt, ' Reset hdmpfmt to 10 as GRIB format'
    hdmpfmt = 10
  END IF

  IF ( lbcopt == 2 ) THEN
    lexbc = 1
  END IF
!
!-----------------------------------------------------------------------
!
!  If we are restart this run from data from a previous run,
!  check if restart data is correct and bring in the data.
!
!  In this case there is no need to initialize the data structure
!
!-----------------------------------------------------------------------
!
  IF (rstart) THEN
!
!-----------------------------------------------------------------------
!
!  Read in restart data at time timestart for all grids.
!
!-----------------------------------------------------------------------
!
    CALL rstrdwr(1,rstime)

    dtbig = possk(1)   ! Use the same big time step in restart file
!
!-----------------------------------------------------------------------
!
!  Initialize arrays that store the lookup table data.
!
!-----------------------------------------------------------------------
!
    CALL initlktb

  END IF

  IF ( lexbc == 1 ) THEN
    CALL setexbcptr( nxc,nyc,nzc )
  END IF

  iorder = 2       ! not used
  tol    = 0.1     ! not used
  cut    = 0.1     ! not used
  relcut = 0.1     ! not used
  cdist  = 0.1     ! not used
  bzone  = 0.1     ! not used

  IF (rstart) THEN

    CALL resett

    CALL stst3
!
!-----------------------------------------------------------------------
!
!  If we're not restarting, but rather starting from the coarse
!  grid only, we need to initialize the data structure and perform
!  the initial gridding.  Note, the timestep sizes and
!  temporal refinement ratio can be changed on the fly.
!  The spatial step sizes and refinement ratios can also be
!  changed on the fly.  Note that here we only change these
!  values in the adaptive grid structure, the user is responsible
!  for making any changes in the real and integer solver constants
!  through the use of resetcon and the reset-constants file
!
!-----------------------------------------------------------------------
!
  ELSE

    CALL stgrid
    CALL stst1
    CALL domain( nxc,nyc,nzc ,tstart )
    CALL stst2

  END IF

  dtbase = possk(1)
!
!-----------------------------------------------------------------------
!
!  set new fine grids if needed
!
!-----------------------------------------------------------------------
!
  IF(grdsrt) THEN

    IF(verbose5)PRINT*,'regrid at restart. To call regrid etc...'

    CALL regrid( levfix,.true. )

  END IF

!
!-----------------------------------------------------------------------
!
!  Set the io times to number of coarse grid timesteps
!
!-----------------------------------------------------------------------
!
  iout(1) = nint(tfmtprt/possk(1))
  IF ( iout(1) == 0 ) iout(1)=-1

  iout(2) = nint(thisdmp/possk(1))
  IF ( iout(2) == 0 ) iout(2)=-1

  iout(3) = nint(trstout/possk(1))
  IF ( iout(3) == 0 ) iout(3)=-1

  iout(4) = nint(tmaxmin/possk(1))
  IF ( iout(4) == 0 ) iout(4)=-1

  iout(5) = nint(tplots/possk(1))
  IF ( iout(5) == 0 ) iout(5)=-1

  iout(6) = nint(tstrtdmp/possk(1))
  IF ( iout(6) == 0 ) iout(6)=-1
!
!-----------------------------------------------------------------------
!
!  Reset certain parameters using values from the solver for the
!  base grid. The values read in from the interface input file
!  are superceded.
!
!-----------------------------------------------------------------------
!
  nx = node(5,mptr)
  ny = node(6,mptr)
  nz = node(14,mptr)

  ii = igtint(mptr,1)
  ir = igtrel(mptr,1)
!
!  Retrieve constant values from the constant arrays, and feed
!  the values into parameters in globcst.inc of ARPS model.
!
  CALL getcnts(nx,ny,nz, a(ii),nsint, a(ir),nsreal)

  IF ( lbcopt /= 2 .AND. lexbc == 1 ) THEN
    WRITE (6,'(a/a)')                                                   &
        ' Warning: model configured NOT to use the external boundary',  &
        ' condition. Some space may be wasted.'
  ELSE IF ( lbcopt == 2 .AND. lexbc /= 1 ) THEN
    WRITE (6,'(a/a/a/a)')                                               &
        ' Error: model configured to use the external boundary ',       &
        ' condition, but no space was allocated to EXBC arrays.',       &
        ' Check parameter lexbc in arpsagri and lbcopt in ',            &
        ' boundary_condition_options. Program stopped.'
    STOP
  END IF

  runname = runnew
  lfnkey  = nmlnthn

  IF ( .NOT.rstart ) THEN

    iout(1) = nfmtprt
    iout(2) = nhisdmp
    iout(3) = nrstout
    iout(4) = nmaxmin
    iout(5) = nplots
    iout(6) = nstrtdmp

  ELSE

    tstart = rstime    ! Use the restart time for tstart
    nfmtprt = iout(1)
    nhisdmp = iout(2)
    nrstout = iout(3)
    nmaxmin = iout(4)
    nplots  = iout(5)
    nstrtdmp= iout(6)
    tfmtprt = MAX(iout(1)*dtbase, 0.0)
    thisdmp = MAX(iout(2)*dtbase, 0.0)
    trstout = MAX(iout(3)*dtbase, 0.0)
    tmaxmin = MAX(iout(4)*dtbase, 0.0)
    tplots  = MAX(iout(5)*dtbase, 0.0)
    tstrtdmp= MAX(iout(6)*dtbase, 0.0)
    tstop = tstopnew
!
! note: here we changed the parameters for the base grid only.
!    The parameters for the other grids should also be
!    changed accordingly. The full implementation of this function
!    is done in the change_constant step.
!
  END IF

  WRITE(6,'(1x,a,i3)')    'For grid ',mptr
  WRITE(6,'(1x,a,i3)')    'nfmtprt: ',nfmtprt
  WRITE(6,'(1x,a,i3)')    'nhisdmp: ',nhisdmp
  WRITE(6,'(1x,a,i3)')    'nrstout: ',nrstout
  WRITE(6,'(1x,a,i3)')    'nmaxmin: ',nmaxmin
  WRITE(6,'(1x,a,i3)')    'nplots:  ',nplots
  WRITE(6,'(1x,a,f10.3)') 'tfmtprt: ',tfmtprt
  WRITE(6,'(1x,a,f10.3)') 'thisdmp: ',thisdmp
  WRITE(6,'(1x,a,f10.3)') 'trstout: ',trstout
  WRITE(6,'(1x,a,f10.3)') 'tmaxmin: ',tmaxmin
  WRITE(6,'(1x,a,f10.3)') 'tplots:  ',tplots
  WRITE(6,'(1x,a,f10.3)') 'tstop:   ',tstop
  WRITE(6,'(1x,a,f10.3)') 'tstrtdmp:',tstrtdmp

  CALL strcnts(nx,ny,nz, a(ii),nsint, a(ir),nsreal)

  nstart= nint( tstart/possk(1) )
  nstop = nint( tstop /possk(1) )
!
!-----------------------------------------------------------------------
!
!  reset the constants in the solver constants data if necessary
!
!-----------------------------------------------------------------------
!
!  if( resetcon ) call chgcst( cstfile )
!
!-----------------------------------------------------------------------
!
!  Print out some I/O control information:
!
!-----------------------------------------------------------------------
!
  WRITE(6,*) ' Adaptive model run ',runname(1:lfnkey)

  IF( rstart ) WRITE(6,*) ' Restart from run ',runold(1:nmlntho)

  WRITE(6,102) tstart,nstart, tstop,nstop

  WRITE(6,103) tol,iorder,kcheck,bzone,cut,relcut,                      &
               cdist,mxnest,intrat,intratt,                             &
               hxposs(1),hyposs(1),possk(1)
!
!-----------------------------------------------------------------------
!
!  output the grid structure
!
!-----------------------------------------------------------------------
!
  CALL outtre( mstart,.true.,.false. )
!
!-----------------------------------------------------------------------
!
!  output initial data
!
!-----------------------------------------------------------------------
!
  CALL usrout1( 1,lfine,1 )
!
!-----------------------------------------------------------------------
!
!  Dump out restart data at initial time.
!
!-----------------------------------------------------------------------
!

  IF(.false.) THEN
    WRITE(6,'(''  DUMPING RESTRART DATA AT INITIAL TIME='',F10.2)')     &
          tstart
    CALL rstrdwr(2,tstart)
  END IF
!
!-----------------------------------------------------------------------
!
!  Call tick to perform time integration ...
!
!-----------------------------------------------------------------------
!

  CALL tick( iout,nstart, nstop )

!
!-----------------------------------------------------------------------
!
!  End of time integration.
!
!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!
!  dump out restart data if wanted
!
!-----------------------------------------------------------------------

  IF (rstdump) THEN
    WRITE(6,'(''  DUMPING RESTART DATA '')')
    timeend = AMAX1(0.,possk(1)*nstop)
    CALL rstrdwr(2,timeend)
  END IF
!
!-----------------------------------------------------------------------
!
!  close graphics stuff
!
!-----------------------------------------------------------------------
!
!  call xgrend
!

!
!-----------------------------------------------------------------------
!
!  Print out some run time information.
!
!-----------------------------------------------------------------------
!
  WRITE(6,*) '   alloc dimension was     ',lstore
  WRITE(6,*) '   maximum alloc used was  ',ihighwater
  WRITE(6,104)
!
!-----------------------------------------------------------------------
!
!  Print the cpu statistics
!
!-----------------------------------------------------------------------
!
  cpu_simulation = f_cputime() - cpu_simulation
  CALL prtcpu(1,6)

  102   FORMAT(//,                                                      &
               ' timing parameters ',//,                                &
               '    start time       ',e12.5,' number of steps ',i6,/,  &
               '    stop time        ',e12.5,' number of steps ',i6,/,  &
               '    output1 interval ',e12.5,' number of steps ',i6,/,  &
               '    output2 interval ',e12.5,' number of steps ',i6,/,  &
               '    output3 interval ',e12.5,' number of steps ',i6,/,  &
               '    output4 interval ',e12.5,' number of steps ',i6,/,  &
               '    output5 interval ',e12.5,' number of steps ',i6 )

  103   FORMAT(///,                                                     &
               ' mesh refinement parameters ',//,                       &
               ' error tol            ',e12.5,/,                        &
               ' order of integrator     ',i9,/,                        &
               ' error checking interval ',i9,/,                        &
               ' buffer zone size     ',e12.5,/,                        &
               ' volume cutoff ratio  ',e12.5,/,                        &
               ' relative vol. cutoff ',e12.5,/,                        &
               ' cluster seperation   ',e12.5,/,                        &
               ' max. refinement level   ',i9,/,                        &
               ' refine ratio, spatial   ',i9,/,                        &
               ' refine ratio, temporal  ',i9,/,                        &
               ' base delta x         ',e12.5,/,                        &
               ' base delta y         ',e12.5,/,                        &
               ' base delta t         ',e12.5,//)
  104   FORMAT(' ------  END OF MR INTEGRATION --------  ')
  105   FORMAT(//,'  changing real and integer constants in the ',/,    &
                  '  solver constants file ',/)

  STOP
END PROGRAM agri
!

SUBROUTINE rstrdwr(i,time)
!
!-----------------------------------------------------------------------
!
!  Reads in or dumps out all the data for a restart
!
!  i = 1 for read in stuff
!    = 2 for write out stuff
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'agrigrid.inc'
  INCLUDE 'manage.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agricst.inc'
  INCLUDE 'globcst.inc'

  REAL :: time
  INTEGER :: i
  CHARACTER (LEN=80) :: filenm,runnami,filenmi,machsti
  REAL :: dum, timei, hzposs
  INTEGER :: igtunit, nmlnthi, iunit ,length
  INTEGER :: nsinti,nsreali,                                            &
          nx1di,ny1di,nz1di,                                            &
          nxy2di,nxz2di,nyz2di,                                         &
          nxyz3di,nexbc3di,lstorei
!
!-----------------------------------------------------------------------
!
!  To read in restart data.
!
!-----------------------------------------------------------------------
!
  IF(i == 1) THEN

    WRITE(6,*) ' restart read from run ',runold(1:nmlntho),             &
               ' at time ',time
!
    CALL mkflnm( runold(1:nmlntho),' ','.t',time,0,filenm,length )
    iunit = igtunit( dum )
    CALL fexist( filenm,length,.true.,.false. )

    OPEN( UNIT=iunit,FILE=filenm(1:length),FORM='unformatted',          &
          STATUS='old')

    IF(verbose3) WRITE(6,*) ' constants file is ',filenm(1:length),     &
                 ', io through unit ',iunit

    REWIND(iunit)
!
! read in the rnuname and
! variable numers and see if they're the same as
! what this program wants
!
    READ(iunit) runnami,filenmi,machsti,nmlnthi,timei
    PRINT*,' runnami,filenmi,machsti,nmlnthi,timei',                    &
             runnami,filenmi,machsti,nmlnthi,timei

    IF(  (runnami(1:nmlnthi) /= runold(1:nmlntho)) .OR.                 &
           (filenmi(1:length)  /= filenm(1:length))  .OR.               &
           (nmlntho /= nmlnthi )                         ) THEN

      WRITE(6,'(/2x,a,/2x,a/)' )                                        &
           'Warning: Inconsistencies found between restart',            &
           'file internal info and the input specifications.'

      WRITE(6,*)  '  runname error in restart '
      WRITE(6,*)  '  input runname           ',runold(1:nmlntho)
      WRITE(6,*)  '  runname in restart file ',runnami(1:nmlnthi)
      WRITE(6,*)  '  input file               ',filenm(1:length)
      WRITE(6,*)  '  input file internal name ',filenmi(1:80)
      WRITE(6,*)  '  time given for restart ',time
      WRITE(6,*)  '  time in restart file   ',timei


!      stop

    ELSE
      WRITE(6,*)  ' time in restart file   ',timei
    END IF

    READ(iunit) nsinti,nsreali,                                         &
                nx1di,ny1di,nz1di,                                      &
                nxy2di,nxz2di,nyz2di,                                   &
                nxyz3di,nexbc3di,lstorei

    IF( (nsint  /= nsinti)  .OR. (nsreal /= nsreal)  .OR.               &
          (nx1d   /= nx1di)   .OR. (ny1d   /= ny1di)   .OR.             &
          (nz1d   /= nz1di)   .OR.                                      &
          (nxy2d  /= nxy2di)  .OR. (nxz2d  /= nxz2di)  .OR.             &
          (nyz2d  /= nyz2di)  .OR.                                      &
          (nxyz3d /= nxyz3di) .OR. (nexbc3d /= nexbc3di) ) THEN

      WRITE(6,'('' ERROR - GRID VARIABLES DIFFER IN NUMBER '')')
      WRITE(6,101)  nsint,nsinti,nsreal,nsreali,                        &
                    nx1d,nx1di,ny1d,ny1di,nz1d,nz1di,                   &
                    nxy2d,nxy2di,nxz2d,nxz2di,nyz2d,nyz2di,             &
                    nxyz3d,nxyz3di,nexbc3d,nexbc3di
      STOP
    END IF
!
!  output storage size information
!
    WRITE(6,'(''  STORE SIZE DATA '',                                   &
    &       /,''    the INPUT storage size was   '',i12,                &
    &       /,''    the PROGRAM storage size is  '',i12 )')             &
           lstorei,lstore
!
!  stop if storage size not large enough
!
    IF(lstorei > lstore) THEN
      WRITE(6,'(''  STORAGE TOO SMALL, ERROR STOP '')')
      STOP
    END IF
!
!  now bring in all the restart data
!  first the nodal data
!
    IF(verbose3) WRITE(6,'(''  READING IN NODAL DATA '')')
    READ(iunit) rnode,node,lstart,newstl,llist,lback,tol,bzone,         &
                mstart,ndfree,intrat,intratt,                           &
                lfine,kcheck,mxnest,hxposs,                             &
                hyposs,hzposs,possk,ncheck,levmlt
!
!  management data
!
    IF(verbose3) WRITE(6,'(''  READING IN MANAGEMENT DATA '')')
    READ(iunit) lfree,lenf,idimf
!
!  read in the grid description data
!
    IF(verbose3) WRITE(6,'(''  READING IN GRID DESCRIPTION DATA '')')

    READ(iunit) stgxy,stgxz,stgyz,stgxyz,stgexbc,                       &
                idmxy,idmxz,idmyz,idmxyz,idmexbc,                       &
                ipkxy,ipkxz,ipkyz,ipkxyz,ipkexbc,                       &
                inixy,inixz,iniyz,inixyz,iniexbc,                       &
                iupxy,iupxz,iupyz,iupxyz,iupexbc,                       &
                ibdxy,ibdxz,ibdyz,ibdxyz,ibdexbc
!
!  read in the grid storage information
!
    IF(verbose3) WRITE(6,'(''  READING IN GRID STORAGE DATA '')')
    READ(iunit) ipint,ipreal,ips1d,                                     &
                ipx,ipy,ipz,ipxy,ipxz,ipyz,ipxyz,ipexbc
    READ(iunit) ntemp

    CLOSE(UNIT=iunit,STATUS='keep')
    IF(verbose3)   WRITE(6,*) ' closed ',filenm(1:length),              &
                     ' on unit ',iunit
    CALL retnunit( iunit )
!
!  finally, read in the grid data
!
    IF(verbose3) WRITE(6,'(''  READING IN GRID DATA '')')
    CALL readall( time )
    WRITE(6,'(''  RESTART READ COMPLETE '')')

  ELSE IF( i == 2 ) THEN

!
!-----------------------------------------------------------------------
!
!  To dump out restart data.
!
!-----------------------------------------------------------------------
!
    WRITE(6,'(1x,a,a,a,f10.2,a)')                                       &
         ' Restart dump for run ',runname(1:lfnkey),                    &
         ' at time ',time,' (s).'

    CALL mkflnm( runname(1:lfnkey),' ','.t',time,0,filenm,length )

    iunit = igtunit( dum )
    CALL fexist( filenm,length,.false.,.true. )

    WRITE(6,'(1x,a)')                                                   &
         ' Restart dump file: '//filenm(1:length)

    OPEN( UNIT=iunit,FILE=filenm(1:length),FORM='unformatted',          &
          STATUS='new')

    IF(verbose3) WRITE(6,*) ' constants file is ',filenm(1:length),     &
                 ', io through unit ',iunit

    IF(verbose3) WRITE(6,'('' DUMPING RUNNAME AND STORE PARAMETERS '')')

!    machst = ' CRAY YMP at NCAR '
    machst = ' UNKNOWN '

    WRITE(iunit) runname,filenm,machst,lfnkey,time
    WRITE(iunit) nsint,nsreal,                                          &
                 nx1d,ny1d,nz1d,                                        &
                 nxy2d,nxz2d,nyz2d,                                     &
                 nxyz3d,nexbc3d,ntemp

    IF(verbose3) WRITE(6,'(''  DUMPING NODEL DATA '')')
    WRITE(iunit) rnode,node,lstart,newstl,llist,lback,tol,bzone,        &
                 mstart,ndfree,intrat,intratt,                          &
                 lfine,kcheck,mxnest,hxposs,                            &
                 hyposs,hzposs,possk,ncheck,levmlt

    IF(verbose3) WRITE(6,'(''  DUMPING MANAGEMENT DATA '')')
    WRITE(iunit) lfree,lenf,idimf

    IF(verbose3) WRITE(6,'(''  DUMPING GRID DESCRIPTION DATA '')')
    WRITE(iunit) stgxy,stgxz,stgyz,stgxyz,stgexbc,                      &
                 idmxy,idmxz,idmyz,idmxyz,idmexbc,                      &
                 ipkxy,ipkxz,ipkyz,ipkxyz,ipkexbc,                      &
                 inixy,inixz,iniyz,inixyz,iniexbc,                      &
                 iupxy,iupxz,iupyz,iupxyz,iupexbc,                      &
                 ibdxy,ibdxz,ibdyz,ibdxyz,ibdexbc

    IF(verbose3) WRITE(6,'(''  DUMPING GRID STORAGE DATA '')')
    WRITE(iunit) ipint,ipreal,ips1d,                                    &
                 ipx,ipy,ipz,ipxy,ipxz,ipyz,ipxyz,ipexbc
    WRITE(iunit) ntemp

    CLOSE(UNIT=iunit,STATUS='keep')
!
!  For PSC Cray, save the restart file in FAR, and
!  delete the file if saved sucessfully.
!

    IF(verbose3) WRITE(6,*) ' closed ',filenm(1:length),                &
                  ' on unit ',iunit
    CALL retnunit( iunit )

    IF(verbose3) WRITE(6,'(''  DUMPING GRID DATA '')')

    CALL dmpall( time )

    WRITE(6,'(''  RESTART DUMP COMPLETE '')')

  END IF
!
!  we're done here
!
  RETURN

  101   FORMAT('            in dataset    in program ',                 &
               '  nsint   ',i10,3X,i10,                                 &
               '  nsreal  ',i10,3X,i10,                                 &
               '  nx1d    ',i10,3X,i10,                                 &
               '  ny1d    ',i10,3X,i10,                                 &
               '  nz1d    ',i10,3X,i10,                                 &
               '  nxy2d   ',i10,3X,i10,                                 &
               '  nxz2d   ',i10,3X,i10,                                 &
               '  nyz2d   ',i10,3X,i10,                                 &
               '  nxyz3d  ',i10,3X,i10,                                 &
               '  nexbc3d ',i10,3X,i10                   )

END SUBROUTINE rstrdwr

SUBROUTINE readst(a,m,iunit)
  REAL :: a(m)
  READ(iunit) a
  RETURN
END SUBROUTINE readst

SUBROUTINE writst(a,m,iunit)
  REAL :: a(m)
  WRITE(iunit) a
  RETURN
END SUBROUTINE writst
