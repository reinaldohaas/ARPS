!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INITASSIM                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initassim(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the recovery/assimilation control parameters. Most of
!  them are read in from an input file. This subroutine also produces
!  a log file which can be used as the input file for replicating a
!  particular recovery/assimilation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Steven Lazarus and Alan Shapiro
!
!  2/20/1993.
!
!  MODIFICATION HISTORY:
!
!  01/03/96 (Limin Zhao)
!  Added control parameters for blending ADAS data with Doppler
!  retrieval(v_err, hfillerr, adas_err).
!
!  11/15/96 (Limin Zhao)
!  Added control parameters for sensitivity tests(bc_opt, itest,
!  ig, ip_wgt, nadas).
!
!-------------------------------------------------------------------------
!
!  INPUT:
!
!    nx        Number of grid points in the x-direction
!    ny        Number of grid points in the y-direction
!    nz        Number of grid points in the z-direction
!
!  OUTPUT:
!
!    Parameters declared in include file 'ASSIM.INC':
!
!    assimnm   The name of this job (a string at least 6 characters
!              long).
!
!    assimdat  File names of the input data.
!
!    dtyp      Input data file type:
!
!              = 0, Model History file format. Contains data at 1 time level.
!              = 1, Processed Lincoln Lab data. (See ARPSPROC.F)
!              = 2, Other. (User specified)
!
!     recovopt Dynamic retrieval option:
!
!              = 0, Do NOT retrieve p' or T'
!              = 1, Retrieve p' and T'
!
!     varopt Variational adjustment option:
!
!              = 0, NO variational adjustment
!              = 1, Perform a variational adjustment
!
!     insrtopt Insertion option:
!
!              = 0, Do NOT insert velocities
!              = 1, Insert velocities
!
!     xshift   x coordinate of lower left hand corner point of
!              model grid with repsect to the radar.
!
!     yshift   y coordinate of lower left hand corner point of
!              model grid with repsect to the radar.
!
!     zshift   z coordinate of lower left hand corner point of
!              model grid with repsect to the radar.
!
!     dirnam  Directory from which the input data is read
!     ldirnm   Length of directory name where the data resides
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations. (Local Variables)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz      ! The number of grid points in 3 directions
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!

  INTEGER :: i

  CHARACTER (LEN=256  ) :: logfn  ! A string used as the log file filename.
  CHARACTER (LEN=128  ) :: stuff  ! A string to 'stuff' in unneeded comment lines.

  INTEGER :: llogfn      ! The length of the log file filename.
  INTEGER :: logfunt     ! FORTRAN unit number for log file output.
  INTEGER :: lenstr      ! Length of a string
  INTEGER :: istat       ! Flag set by open statement on the status
                         ! of file opening
  INTEGER :: inparm      ! I/O unit for the ASSIM.INPUT file
  INTEGER :: outparm     ! I/O unit for writing assimilation parameters

  LOGICAL :: iexist      ! Flag set by inquire statement for file
                         ! existence
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'assim.inc'
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  namelist Declarations:
!
!-----------------------------------------------------------------------
!
  NAMELIST /assimilation/ assimopt, assimnm,dirnam,dtyp,                &
                          recovopt,varopt,insrtopt,                     &
                          xshift,yshift,zshift,                         &
                          ip_wgt,v_err,vfill_err,adas_err,spval,        &
                          nvf,assimdat,nadas,adasdat,                   &
                          itest,bc_opt,ig,                              &
                          voltim1,voltim2,voltim3,                      &
                          radid,latrad,lonrad,elvrad
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
  EXTERNAL assimpar
  EXTERNAL strlnth
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL getunit(outparm)

!  open(unit=outparm,file='assimbat.out',form='formatted',
!    :     status='unknown',iostat=istat)

!  IF( istat.ne.0) THEN

!    write(outparm,'(/1x,a,i2,/1x,a/)')
!    :    'Error occured when opening output file ASSIMBAT.OUT
!    :     in FORTRAN unit ',outparm,'Program stopped in ASSIMBAT.'
!    STOP

!  ENDIF

  outparm = 6

  WRITE(outparm,'(/ 15(/5x,a)//)')                                       &
     '###############################################################',  &
      '###############################################################', &
      '#####                                                     #####', &
      '#####                  Welcome to                         #####', &
      '#####                                                     #####', &
      '#####      The recovery/assimilation input parameter      #####', &
      '#####               suboutine ASSIMBAT.F                  #####', &
      '#####                                                     #####', &
      '#####                                                     #####', &
      '#####                      Developed by                   #####', &
      '#####      Center for Analysis and Prediction of Storms   #####', &
      '#####                 University of Oklahoma              #####', &
      '#####                                                     #####', &
      '###############################################################', &
      '###############################################################'
!
!
!-----------------------------------------------------------------------
!
!  Read in the values of the recovery/assimilation parameters:
!
!-----------------------------------------------------------------------
!

!  CALL getunit(inparm)

!  open(unit=inparm,file='assim.input',form='formatted',
!   open(unit=inparm,file='assim.input.2',form='formatted',
!    :     status='old',iostat=istat)

!  IF( istat.ne.0) THEN

!    write(outparm,'(/1x,a,i2,/1x,a/)')
!    :  'Error occured when opening initial input file ASSIM.INPUT
!    :    using FORTRAN unit ',inparm,'Program stopped in ASSIMBAT.'
!    STOP

!  ENDIF

  inparm = 5

  WRITE(outparm,'(4(/5x,a),/6(/5x,a)/)')                                &
      'The recovery/assimilation is run concurrently with the ARPS',    &
      'model by inputing a number of control parameters from an',       &
      'external file, assim.input. The parameters are input via the',   &
      'screen for interactive jobs.',                                   &
      'After completing the parameter input, a log file is produced,',  &
      'which can then be used as the input file to replicate the',      &
      'job. This file is named assimnm.log.nn, where assimnm is a',     &
      'prefix for the output files that are associated with the',       &
      'assimilation and nn is a number appended to the file name when', &
      'file assim.log.nn-1 already exists.'


!
!-----------------------------------------------------------------------
!
!  Read in comment lines for this job and the name of assimilation log
!  file - designated by a string at least 6 characters long.
!
!-----------------------------------------------------------------------
!
  WRITE(outparm,'(/2x,a,/3(/5x,a),a)')                                  &
            'COMMENT LINES AND JOB NAME',                               &
            'You can provide the model with a few lines of comments',   &
            'which designate the nature of this model run. You need',   &
            'to tell the model first how many lines of comments you',   &
            'have.'

  WRITE(outparm,'(/5x,a/)')                                             &
       'Input the number of comment lines (0 to 50):'

  READ (inparm,assimilation)

  WRITE(outparm,'(/5x,a,/5x,a,/5x,a,/5x,a/)')                           &
       'The assimilation requires a 6 character long string to',        &
       'construct the names of output files. This string is the',       &
       'the name of the log file. Please input such a string',          &
       '(at least 6 characters long):'

  WRITE(outparm,'(/5x,a,i3)')                                           &
      'The ARPS data assimilation option is: ', assimopt

  WRITE(outparm,'(/5x,a,a)')                                            &
      'The name of this run is: ', assimnm

  IF( assimnm == runname) THEN
    WRITE (outparm,'(3(/5x,a))')                                        &
        'WARNING: The assimilation name matches that of the run name',  &
        'Job aborted in subroutine INITASSIM.'
    STOP

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in the number of input files (minimum of 3).
!
!-----------------------------------------------------------------------
!
  WRITE(outparm,'(/3(/5x,a))')                                          &
       'In order to perform an assimilation, there must',               &
       'be a minimum of 3 input files at 3 different times.',           &
       'Please set the number of input data files >=3'

  WRITE(outparm,'(5x,a,i4)') 'Input was ',nvf
  IF ( nvf < 3 ) THEN
    WRITE (outparm,'(/5x,a,i3,a))')                                     &
        'WARNING: The number of input velocity files,', nvf,            &
        ' is less than 3.'
    WRITE (outparm,'(/5x,a))')                                          &
        'Job aborted in subroutine INITASSIM.'
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in the input data file type
!
!  dtyp= 0, input files are model history data files
!          = 1, input files are processed Lincoln Lab data files
!          = 2, other (user specified)
!
!  For all options, the names of input files need to be provided.
!
!-----------------------------------------------------------------------
!
  WRITE(outparm,'(/2x,a,//5x,a,//5x,a,/5x,a,/5x,a,//5x,a)')             &
            'INPUT DATA FILE TYPE:',                                    &
            'Two options:',                                             &
            '  0  input files are model history data files ',           &
            '  1  input files are processed Lincoln Lab data files ',   &
            '  2  input files are other (user specified) ',             &
            'Select input data file type (see above):'

  WRITE(outparm,'(5x,a,i4)') 'input was ',dtyp

  IF( dtyp /= 0.AND.dtyp /= 1.AND.dtyp /= 2) THEN
    WRITE (outparm,'(3(/5x,a))')                                        &
        'dtyp must equal either 0,1 or 2',                              &
        'Please set dtyp to 0, 1, or 2 prior to running the model.',    &
        'Job aborted in subroutine ASSIMBAT.'
    STOP

  END IF
!
!
!-----------------------------------------------------------------------
!
!  Dynamic retreival flag
!
!  recovopt = 0, Do NOT retrieve p, T
!           = 1, Retrieve p, T
!
!-----------------------------------------------------------------------
!
  WRITE(outparm,'(/2x,a,//5x,a,//5x,a,/5x,a,//5x,a)')                   &
            'DYNAMIC RETRIEVAL FLAG:',                                  &
            'Two options:',                                             &
            '  0  no insertion, no retrieval',                          &
            '  1  insert velocities/retrieve p, T',                     &
            'Select retrieval option (see above):'
!
!  read (inparm,*) recovopt
  WRITE(outparm,'(/5x,a,i4/)') 'Input was ', recovopt

  IF( recovopt /= 0.AND.recovopt /= 1) THEN
    WRITE (outparm,'(3(/5x,a))')                                        &
        'recovopt must equal either 0,or 1',                            &
        'Please set recovopt to 0,or 1 prior to running the model.',    &
        'Job aborted in subroutine ASSIMBAT.'
    STOP

  END IF

!  read(inparm,*) stuff
!  read(inparm,*) stuff
!
!-----------------------------------------------------------------------
!
!  Variational adjustment flag
!
!  varopt  = 0, NO variational adjusment
!             = 1, Perform variational adjustment
!
!-----------------------------------------------------------------------
!
  WRITE(outparm,'(/2x,a,//5x,a,//5x,a,/5x,a,//5x,a)')                   &
            'VARIATIONAL ADJUSTMENT FLAG:',                             &
            'Two options:',                                             &
            '  0  NO variational adjustment ',                          &
            '  1  Perform variational adjustment',                      &
            'Select variational adjustment option (see above):'
!
!  read (inparm,*) varopt
  WRITE(outparm,'(/5x,a,i4/)') 'Input was ', varopt

  IF( varopt /= 0.AND.varopt /= 1) THEN
    WRITE (outparm,'(3(/5x,a))')                                        &
        'varopt must equal either 0,or 1',                              &
        'Please set varopt to 0,or 1 prior to running the model.',      &
        'Job aborted in subroutine ASSIMBAT.'
    STOP

  END IF

!  read(inparm,*) stuff
!  read(inparm,*) stuff
!
!-----------------------------------------------------------------------
!
!  Insertion flag
!
!  insrtopt = 0, NO insertion
!           = 1, Perform an insertion
!
!-----------------------------------------------------------------------
!
  WRITE(outparm,'(/2x,a,//5x,a,//5x,a,/5x,a,//5x,a)')                   &
            'INSERTION FLAG:',                                          &
            'Two options:',                                             &
            '  0  NO insertion',                                        &
            '  1  Perform an insertion',                                &
            'Select insertion option (see above):'
!
!  read (inparm,*) insrtopt
  WRITE(outparm,'(/5x,a,i4/)') 'Input was ', insrtopt

  IF( insrtopt /= 0.AND.insrtopt /= 1) THEN
    WRITE (outparm,'(3(/5x,a))')                                        &
        'insrtopt must equal either 0,or 1',                            &
        'Please set insrtopt to 0,or 1 prior to running the model.',    &
        'Job aborted in subroutine ASSIMBAT.'
    STOP

  END IF

!
!-----------------------------------------------------------------------
!
!  If the recovery/assimilation options are turned off, exit.
!
!-----------------------------------------------------------------------
!

  IF (recovopt == 0.AND.varopt == 0.AND.insrtopt == 0) RETURN

!
!-----------------------------------------------------------------------
!
!  Input the position of the lower left hand corner of the cartesian
!  model grid with respect to the radar (xshift,yshift,zshift).
!
!-----------------------------------------------------------------------
!
!  read(inparm,*) stuff
!  read(inparm,*) stuff

  WRITE(outparm,'(2(/5x,a))')                                           &
         'Input the x-ccordinate of the lower left hand corner of:',    &
         'Cartesian model grid:'

!  read (inparm,*) xshift
  WRITE(outparm,'(/5x,a,f8.0/)') 'Input was ', xshift

  WRITE(outparm,'(2(/5x,a))')                                           &
         'Input the y-ccordinate of the lower left hand corner of:',    &
         'Cartesian model grid:'

!  read (inparm,*) yshift
  WRITE(outparm,'(/5x,a,f8.0/)') 'Input was ', yshift

  WRITE(outparm,'(2(/5x,a))')                                           &
         'Input the z-ccordinate of the lower left hand corner of:',    &
         'Cartesian model grid:'

!  read (inparm,*) zshift
  WRITE(outparm,'(/5x,a,f8.0/)') 'Input was ', zshift
!
!-----------------------------------------------------------------------
!
!    Blending control parametrs.
!
!-----------------------------------------------------------------------
!
!  read (inparm,*) stuff
!
  WRITE(outparm,'(1(/5x,a))')                                           &
         'Input the option for determine the blending weights:'

!  read (inparm,*) ip_wgt
  WRITE(outparm,'(/5x,a,i4/)') 'Input was ', ip_wgt
!
  WRITE(outparm,'(1(/5x,a))')                                           &
         'Input the variance error for the retrival velocity:'

!  read (inparm,*) v_err
  WRITE(outparm,'(/5x,a,f8.0/)') 'Input was ', v_err

  WRITE(outparm,'(1(/5x,a))')                                           &
      'Input the variance error for the hole-fill velocity:'

!  read (inparm,*) vfill_err
  WRITE(outparm,'(/5x,a,f8.0/)') 'Input was ', vfill_err

  WRITE(outparm,'(1(/5x,a))')                                           &
         'Input the variance error for adas background:'

!  read (inparm,*) adas_err
  WRITE(outparm,'(/5x,a,f8.0/)') 'Input was ', adas_err

  WRITE(outparm,'(1(/5x,a))')                                           &
      'Input the bad data flag:'

!  read (inparm,*) spval
  WRITE(outparm,'(/5x,a,f8.0/)') 'Input was ',spval
!
!-----------------------------------------------------------------------
!
!  Assimilation input parameter:
!
!  Input the name of the directory from which the data
!  files will be read:
!
!-----------------------------------------------------------------------
!

  WRITE(outparm,'(3(/5x,a),/5x,a,a,/5x,a/)')                            &
      'You can redirect the input of data files from a separate',       &
      'directory rather from than the current work directory by',       &
      'specifying the full path name of such a directory.',             &
      'If you give a blank string the current work directory ',         &
      'is assumed.',                                                    &
      'The directory for input is:'

!  read (inparm,*) dirnam

  ldirnm = 256
  CALL strlnth( dirnam, ldirnm)

  IF( ldirnm == 0 ) THEN
    dirnam = '.'
    ldirnm = 1
  END IF

  IF( dirnam(1:ldirnm) /= ' ') THEN
!
!-----------------------------------------------------------------------
!
!  Check if the specified input directory exists, if not,
!  abort the job.
!
!-----------------------------------------------------------------------
!
    INQUIRE(FILE=dirnam(1:ldirnm),EXIST=iexist)

    IF( .NOT.iexist ) THEN

      WRITE(outparm,'(/5x,a,2(/5x,a))')                                 &
          'Specified input directory '//dirnam(1:ldirnm)//              &
          ' not found.',                                                &
          'Please create it before starting the model.',                &
          'Job aborted in subroutine ASSIMBAT.'
      STOP

    END IF

    WRITE(outparm,'(/5x,a,a)')                                          &
        'Input files will be in directory ',                            &
        dirnam(1:ldirnm)//'.'

  ELSE

    WRITE(outparm,'(/5x,a)')                                            &
        'Data files are in the current work directory.'

  END IF

!  read(inparm,*) stuff
!
!-----------------------------------------------------------------------
!
!    Input grid/base state file name (for model data ingest only)
!
!-----------------------------------------------------------------------
!
!  write(outparm,'(/a/a/a)')
!    :'Please give the name of the file containing the grid and base',
!    :'state array data. This file will not be read if the grid and ',
!    :'base state arrays are present in each time-dependent data file.'

!  read(inparm,*) gbfile

!  lenstr = 256
!  CALL strlnth( gbfile, lenstr)
!  write(outparm,'(/5x,a,a)')
!    :'The grid/base state data file to be read in is',
!    : gbfile(1:lenstr)
!
!-----------------------------------------------------------------------
!
!    Input data files.
!
!-----------------------------------------------------------------------
!
  WRITE(outparm,'(//2x,a,//5x,a,/5x,a/)')                               &
       'INITIAL DATA FILES:',                                           &
       'Input the names of the data files used by the',                 &
       'recovery/assimilation (character string in quotes):'
!
  DO i = 1, nvf
!      read(inparm,*) assimdat(i)

    lenstr = LEN(assimdat(i))
    CALL strlnth( assimdat(i), lenstr)
    WRITE(outparm,'(/5x,a,a)')                                          &
        'The data file to be read in is ',                              &
        assimdat(i)(1:lenstr)

  END DO
!
!-----------------------------------------------------------------------
!
!  Read in ADAS background fields.
!
!-----------------------------------------------------------------------
!
  WRITE(outparm,'(1(/5x,a))')                                           &
      'Input the number of background fields:'

!  read (inparm,*) nadas
  WRITE(outparm,'(/5x,a,i4/)') 'Input was ', nadas


  WRITE(outparm,'(//2x,a,//5x,a,/5x,a/)')                               &
       'BACKGROUND DATA FILES:',                                        &
       'Input the names of the data files used by the',                 &
       'recovery/assimilation (character string in quotes):'
!
  DO i = 1, nadas

!      read(inparm,*) adasdat(i)

    lenstr = LEN(adasdat(i))
    CALL strlnth( adasdat(i), lenstr)
    WRITE(outparm,'(/5x,a,a)')                                          &
        'The data file to be read in is ',                              &
        adasdat(i)(1:lenstr)

  END DO
!
!-----------------------------------------------------------------------
!
!  Several options for doing sensitive tests on pressure B.C.,
!  different strategy to do hole-filling, estimation of
!  the perturbation pressure at the lowest levels.
!
!-----------------------------------------------------------------------
!
  WRITE(outparm,'(1(/5x,a))')                                           &
      'Input the option for pressure B.C. in psolver'

!  read (inparm,*) bc_opt
  WRITE(outparm,'(5x,a,i4)') 'Input was ',bc_opt

  WRITE(outparm,'(1(/5x,a))')                                           &
      'Input the option for hole-filling test'

!  read (inparm,*) itest
  WRITE(outparm,'(5x,a,i4)') 'Input was ',itest

  WRITE(outparm,'(1(/5x,a))')                                           &
      'Input the option for estimate the pressure at the lowest level'

!  read (inparm,*) ig
  WRITE(outparm,'(5x,a,i4)') 'Input was ', ig

!
!-----------------------------------------------------------------------
!
!  Print out the input parameters.
!  Write out a log file of model parameters which can be used as
!  the input file to re-run the model.
!
!-----------------------------------------------------------------------
!
  CALL assimpar( nx, ny, nz)

  CALL retunit(outparm)
  CALL retunit(inparm)

  RETURN
END SUBROUTINE initassim

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASSIMPAR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE assimpar(nx,ny,nz)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Print out the input parameters.
!  Write out a log file of model parameters which can be used as
!  the input file to re-run the model.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Steven Lazarus and Alan Shapiro
!
!  2/23/1993.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction
!    ny       Number of grid points in the y-direction
!    nz       Number of grid points in the z-direction
!
!  OUTPUT:
!
!    Parameters declared in include file 'assim.inc':
!
!    assimnm   The name of this job (a string at least 6
!              characters long).
!
!    assimdat  File names of input data.
!
!    nvf  Number of input data files (minimum of 3)
!
!    dtyp  Input data file type:
!
!              = 0, Model History file format. Contains data at 1 time level.
!              = 1, Processed Lincoln Lab data. (See ARPSPROC.F)
!              = 2, Other. (User specified)
!
!     recovopt Dynamic retrieval option:
!
!              = 0, Do NOT retrieve p' or T'
!              = 1, Retrieve p' and T'
!
!     varopt   Variational adjustment option:
!
!              = 0, NO variational adjustment
!              = 1, Perform a variational adjustment
!
!     insrtopt Insertion option:
!
!              = 0, Do NOT insert velocities
!              = 1, Insert velocities
!
!     xshift   x coordinate of lower left hand corner point of
!              model grid with repsect to the radar.
!
!     yshift   y coordinate of lower left hand corner point of
!              model grid with repsect to the radar.
!
!     zshift   z coordinate of lower left hand corner point of
!              model grid with repsect to the radar.
!
!     dirnam  Directory name where input data resides
!     ldirnm  Length of directory name, dirnam
!
!     Parameters declared in include file 'globcst.inc':
!
!     mgrid    Grid identifier
!     nestgrd  Grid nesting flag
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations. (Local Variables)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz      ! The number of grid points in 3 directions
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!

  CHARACTER (LEN=256  ) :: logfn  ! A string used as the log file filename.
  INTEGER :: llogfn      ! The length of the log file filename.
  INTEGER :: logfunt     ! FORTRAN unit number for log file output.
  INTEGER :: lenstr      ! Length of a string
  INTEGER :: lengbf      ! Length of a string
  INTEGER :: istat       ! Flag set by open statement on the status
                         ! of file opening
  INTEGER :: iout        ! I/O channel output of misc. log file data
  LOGICAL :: iexist      ! Flag set by inquire statement for file
                         ! existence
  INTEGER :: i
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'assim.inc'
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
  EXTERNAL strlnth
  EXTERNAL gtlogfn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL getunit(iout)

!  open(unit=iout,file='assimbat.misc',form='formatted',
!    :     status='unknown',iostat=istat)

!  IF( istat.ne.0) THEN

!    write(iout,'(/1x,a,i2,/1x,a/)')
!    :    'Error occured when opening misc. output file ASSIM.MISC
!    :      using FORTRAN unit ',iout,'Program stopped in ASSIMPAR.'
!    STOP

!  ENDIF

  iout = 6

  WRITE(iout,'(/5x,a/)')                                                &
       'Comment lines for this run:'

  WRITE(iout,*) ' '
  WRITE(iout,102) 'nvf =',nvf
  WRITE(iout,102) 'dtyp =',dtyp
  WRITE(iout,102) 'recovopt =',recovopt
  WRITE(iout,102) 'varopt=',varopt
  WRITE(iout,102) 'insrtopt =',insrtopt
  WRITE(iout,104) 'xshift   =',xshift
  WRITE(iout,104) 'yshift   =',yshift
  WRITE(iout,104) 'zshift   =',zshift
  WRITE(iout,*) ' '
  WRITE(iout,102) 'ip_wgt =',ip_wgt
  WRITE(iout,104) 'v_err    =',v_err
  WRITE(iout,104) 'vfill_err=',vfill_err
  WRITE(iout,104) 'adas_err  =',adas_err
  WRITE(iout,104) 'spval    =',spval
  WRITE(iout,*) ' '
  WRITE(iout,102) 'bc_opt    =',bc_opt
  WRITE(iout,102) 'itest =',itest
  WRITE(iout,102) 'ig   =',ig

  101   FORMAT(5X,a,f12.3)
  102   FORMAT(5X,a,i4)
  103   FORMAT(5X,a,e13.3)
  104   FORMAT(5X,a,f8.0)

  WRITE(iout,*) ' '

!
!-----------------------------------------------------------------------
!
!  Now let's write out a log file of this job which can be used
!  as the input file for replicating this run.
!
!  First get a name for the log file:
!
!-----------------------------------------------------------------------
!
  lenstr = LEN(assimnm)
  CALL strlnth( assimnm, lenstr )
  CALL gtlogfn(assimnm(1:lenstr), mgrid, nestgrd, logfn, llogfn)

  CALL getunit(logfunt)

  OPEN (UNIT=logfunt, FILE=logfn(1:llogfn),STATUS='new',                &
              IOSTAT=istat)

  IF(istat /= 0) THEN

    WRITE(iout,'(/5x,a)') 'Error in opening log file ', logfn
    WRITE(iout,'(5x,a/)') 'Job stopped in subroutine ASSIMPAR.'
    STOP

  ELSE

    WRITE(iout,'(/5x,a,i3,a)')                                          &
         'FORTRAN unit ',logfunt,                                       &
         ' will be used for log file output.'

  END IF

  WRITE (logfunt, '(/1x,a)')       '&assimilation'

  WRITE (logfunt, '(3x,a,i4,a)')     'assimopt = ', assimopt, ','

  WRITE(logfunt,'(5x,a,a,a)')                                           &
                       'assimnm = ''',assimnm(1:lenstr),''','

  WRITE(logfunt,'(5x,a,a,a)')                                           &
                     'dirnam = ''',dirnam(1:ldirnm),''','

  WRITE (logfunt, '(5x,a,i4,a)')    'dtyp  = ', dtyp, ','
  WRITE (logfunt, '(5x,a,i4,a)')    'recovopt  = ', recovopt, ','
  WRITE (logfunt, '(5x,a,i4,a)')    'varopt = ', varopt,','
  WRITE (logfunt, '(5x,a,i4,a)')    'insrtopt  = ', insrtopt, ','

  WRITE (logfunt, '(3x,a,f16.4,a)') 'xshift    = ',xshift,    ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'yshift    = ',yshift,    ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'zshift    = ',zshift,    ','

  WRITE (logfunt, '(5x,a,i4,a)')    'ip_wgt  = ',ip_wgt,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'v_err  = ',v_err,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'vfill_err = ',vfill_err, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'adas_err   = ',adas_err,   ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'spval     = ',spval,     ','

  WRITE (logfunt, '(5x,a,i4,a)')    'nvf  = ',nvf,  ','

  DO i=1,nvf
    lenstr = LEN( assimdat(i) )
    CALL strlnth (assimdat(i), lenstr )
    WRITE(logfunt,'(7x,a,i2.2,a,a,a)')                                  &
              'assimdat( ',i,') = ''',assimdat(i)(1:lenstr),''','
  END DO

  WRITE (logfunt, '(5x,a,i4,a)')    'nadas = ',nadas, ','

  DO i=1,nadas
    lenstr = LEN( adasdat(i) )
    CALL strlnth (adasdat(i), lenstr )
    WRITE(logfunt,'(7x,a,i2.2,a,a)')                                    &
               'adasdat( ',i,')  = ''',adasdat(i)(1:lenstr),''','
  END DO

  WRITE (logfunt, '(5x,a,i4,a)')    'itest = ',itest, ','
  WRITE (logfunt, '(5x,a,i4,a)')    'bc_opt    = ',bc_opt,    ','
  WRITE (logfunt, '(5x,a,i4,a)')    'ig   = ',ig,   ','

  WRITE (logfunt, '(3x,a,f16.4,a)') 'voltim1   = ',voltim1,   ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'voltim2   = ',voltim2,   ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'voltim3   = ',voltim3,   ','

  lenstr = LEN( radid )
  CALL strlnth( radid, lenstr )
  WRITE(logfunt,'(5x,a,a,a)')                                           &
                              'radid  = ''',radid(1:lenstr),''','

  WRITE (logfunt, '(3x,a,f16.4,a)') 'latrad    = ',latrad,    ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'lonrad    = ',lonrad,    ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'elvrad    = ',elvrad,    ','

  WRITE (logfunt, '(1x,a)')       '&end'

  CLOSE ( logfunt )
  CALL retunit(logfunt)

  WRITE(iout,'(/5x,a,a,a/)')                                            &
       'Log file ',logfn(1:llogfn),' was produced for this job.'

!  CALL retunit(iout)

  RETURN
END SUBROUTINE assimpar
