SUBROUTINE chgcst( filein )
  CHARACTER (LEN=*) :: filein
  LOGICAL :: existf
  CHARACTER (LEN=1) :: dtype
!-----------------------------------------------------------------
!
!  this routine changes the real and integer constants
!  for the grids.  the input file 'filein' should have
!  the following format for changes:
!
!  each line specifies a change, and each line has
!
!  --->  character, integer1, integer2, value
!
!  character: r or i, for changing a real or an integer constant
!          NOTE:  must be surrounded by quotes, i.e. 'r' or 'i'
!  integer1 : the number of the real or integer const. to change
!  integer2 : the grid level or grid number of the change.
!          positive means all grids on level integer2,
!          negative means change grid abs(integer2),
!          if zero, change all values on all grid - all levels.
!  value    : the new value to assign to the constant
!
!  the constants are defined in agrigrid.inc
!
!  A logical is read in the main program
!  from the unit 7 input file to determine
!  if changes should be made.  If this is true this routine is called.
!  if the changes file is not found or if there is no data in the
!  changes file the program error exits.
!
!-----------------------------------------------------------------

!
!  check the filename
!

  i = INDEX( filein,' ' ) - 1
  IF(i == -1) i = LEN( filein )
  IF(i == 0) THEN
    WRITE(6,*) ' error in filename for changes file ',                  &
               ' length is zero '
    STOP
  END IF
!
  INQUIRE( FILE=filein(1:i),EXIST=existf )
  IF( .NOT. existf ) THEN
    WRITE(6,*) '  constant changes file ',filein(1:i),                  &
               ' does not exist, error exit '
    STOP
  END IF
!
  iunit = igtunit( dum )
  OPEN( UNIT=iunit,FILE=filein(1:i),FORM='formatted',                   &
        STATUS='old'                                    )
!
  WRITE(6,*) ' opened changes file: ',filein(1:i)
  ichanges = 0
!
!  read each change and make it
!

  5    CONTINUE
  ichanges = ichanges + 1
  READ( iunit,*,END=15 ) dtype
  GO TO 16
!
!  if this is the first change and the file is empty,
!  error exit, else exit normally
!
  15   IF(ichanges < 2) THEN
    WRITE(6,*) '  changes file is empty, error exit '
    STOP
  ELSE
    ichanges = ichanges-1
    WRITE(6,*) '  finished making changes, found ',ichanges,            &
               '  in the file ',filein(1:i)
    GO TO 999
  END IF
!
!  make a change
!
  16   CONTINUE
  BACKSPACE( iunit )
  IF( dtype == 'r') THEN
    READ( iunit,* ) dtype, numb, level, value
    CALL chngv( dtype,numb,level,value )
  ELSE IF( dtype == 'i') THEN
    READ( iunit,* ) dtype, numb, level, ivalue
    CALL chngv( dtype,numb,level,ivalue )
  ELSE
    WRITE(6,*) ' error in change file, dtype = ',dtype
    STOP
  END IF
!
  GO TO 5
!
  999  CONTINUE
  CLOSE(UNIT=iunit,STATUS='keep')
  CALL retnunit( iunit )
  RETURN
END SUBROUTINE chgcst
!
!

SUBROUTINE chngv( dtype,numb,level,value )
  INTEGER :: numb,level,gridno
  CHARACTER (LEN=1) :: dtype
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
!
!  this routine changes the real or int const numb to 'value'
!  first check on some error conditions
!
  IF( dtype == 'r' ) THEN
    IF( numb > nsreal ) THEN
      WRITE(6,*) '  numb gt nsreal in chngv.  numb, nsreal are ',       &
                 numb,nsreal
      WRITE(6,*) '  error exit  '
      STOP
    END IF
  ELSE IF( dtype == 'i' ) THEN
    IF( numb > nsint ) THEN
      WRITE(6,*) '  numb gt nsint in chngv.  numb, nsint are ',         &
                 numb,nsint
      WRITE(6,*) '  error exit  '
      STOP
    END IF
  ELSE
    WRITE(6,*) '  error in chngv, dtype is ',dtype
    STOP
  END IF
!
  IF( level > mxnest ) THEN
    WRITE(6,*) '  level is greater than maximum possible level '
    WRITE(6,*) '  in chngrl.  level, mxnest are ',level,mxnest
    WRITE(6,*) '  error exit  '
    STOP
  END IF
!
!  do the change
!
  levst = level
  levend = level
  gridno = 0
  IF( level == 0) THEN
    levst = 1
    levend = lfine
  ELSE IF( level < 0 ) THEN
    gridno = ABS( level )
    levst = 1
    levend = lfine
  END IF
  levelc = levst
!
!  this is the standard loop over the grids in level
!
!  first, for each level
!
  10   IF(levelc > levend) GO TO 999
!
!  now, for each grid
!
  mptr = lstart(levelc)
  5    IF(mptr == 0) GO TO 15
!
  IF ( (mptr == gridno) .OR. (gridno == 0) ) THEN
    IF( dtype == 'r' ) THEN
      in = igtrel( mptr, 1 )
      CALL chgr( mptr,levelc,numb,a(in),value )
    ELSE
      in = igtint( mptr, 1 )
      CALL chgi( mptr,levelc,numb,a(in),value )
    END IF
  END IF
!
  mptr = node(10,mptr)
  GO TO 5
!
! change level when we've exhausted the grids on levelc
!
  15   levelc = levelc + 1
  GO TO 10
!
! we're at the end
!
  999  CONTINUE
  RETURN
END SUBROUTINE chngv
!
!

SUBROUTINE chgr( mptr,level,numb,a,value )
  INTEGER :: numb,level,mptr
  REAL :: value, a(numb)
!
!  make the change
!
  WRITE(6,10) numb,mptr,level,a(numb),value
  a(numb) = value
  10   FORMAT('  change real ',i3,', grid ',i3,', level ',              &
              i2,', from ',e10.3,' to ',e10.3)
  RETURN
END SUBROUTINE chgr
!
!

SUBROUTINE chgi( mptr,level,numb,ia,ivalue )
  INTEGER :: numb,level,mptr,ivalue,ia(numb)
!
!  make the change
!
  WRITE(6,10) numb,mptr,level,ia(numb),ivalue
  ia(numb) = ivalue
  10   FORMAT('  change intg ',i3,', grid ',i3,', level ',              &
              i2,', from ',i10,' to ',i10)
  RETURN
END SUBROUTINE chgi


