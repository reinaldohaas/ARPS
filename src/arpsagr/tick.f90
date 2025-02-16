SUBROUTINE tick( iout,nstart, nstop )

  INCLUDE 'agrialloc.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'agricst.inc'

  DIMENSION intcnt(10),icheck(10),iout(6)
  LOGICAL :: tprint,aprint
  DATA      tprint /.true./
  DATA      aprint/.false./
!
!-----------------------------------------------------------------------
!
!  parameters:
!
!  nstart  = # of coarse grid time steps corresponding to start time.
!  nstop   = # of coarse grid time steps corresponding to stop time.
!  iout    = output interval every 'iout' coarse time steps
!            there are three io options at present
!
!  main driver routine.  controls:
!
!     integration  of all grids.
!     error estimation / regridding
!     output counting
!     updating of fine to coarse grids
!
!  integration strategy is to advance a fine grid until one step of
!  the coarser grid would catch up to it.
!
!  array intcnt counts down from levmlt for the number of times
!  a grid is integrated, to determine who should be integrated next.
!  array icheck also counts (down from ncheck)  to keep track of
!  when that level should have its error estimated/regridding.
!
!-----------------------------------------------------------------------
!
  time      = rnode(20,mstart)
  DO i   = 1, mxnest
    icheck(i) = ncheck(i)
  END DO
!
!-----------------------------------------------------------------------
!
!  start of coarse grid integration loop
!
!-----------------------------------------------------------------------
!
  delt0=possk(mstart)

  DO nstepc=nstart+1,nstop
!
    WRITE(6,'(''  COARSE GRID STEP '',I5)')  nstepc
!
    DO i=1, mxnest
      intcnt(i) = levmlt(i)
    END DO
!
    level     = 1
!
    30     IF(intcnt(level) > 0) GO TO 60
!
!-----------------------------------------------------------------------
!
!    find next level to be integrated
!
!-----------------------------------------------------------------------
!
    40     IF (level <= 1) GO TO 45
    level = level - 1
    IF (intcnt(level) <= 0) GO TO 40
    GO TO 50

    45     IF (level == 1) GO TO 110

    50     lp = level+1

    DO i = lp,lfine
      intcnt(i) = levmlt(i)
    END DO
!
    IF(verbose6) WRITE(6,'(''CALLING UPDATE FOR LEVEL '',I5)') level

    CALL update( level )

    WRITE(6,'(''UPDATE CALLED FOR LEVEL '',I5)') level
!
!-----------------------------------------------------------------------
!
!    regridding  time?
!
!-----------------------------------------------------------------------
!
    60     IF (icheck(level) > 0) GO TO 90
!
!-----------------------------------------------------------------------
!
!    regrid level 'level+1' up to finest level.
!    level 'lbase' stays fixed.
!
!    Not implemented yet.
!
!-----------------------------------------------------------------------
!
    WRITE(6,*) ' calling regrid in tick '
    WRITE(6,*) ' regrid not implemented yet for call here '
    STOP
!
!    lbase = level
!    if (tprint) write(6,101) lbase
!101    format(8h  level ,i5,32h  stays fixed during regridding )
!
!    call regrd1( lbase,.true. )
!
!    call outtre( mstart,.false. )
!
!    maybe finest level in existence has changed. reset counters.
!
!    do 80  i = lbase, lfine
!      icheck(i) = ncheck(i)
!80     continue
!
!-----------------------------------------------------------------------
!
!    done regridding --------------------
!
!-----------------------------------------------------------------------
!

    90     CONTINUE
!
!-----------------------------------------------------------------------
!
!  Integrate all grids at level 'level'.
!
!-----------------------------------------------------------------------
!

    CALL advanc( level)
!
!-----------------------------------------------------------------------
!
!  Done with a level of integration. update counts, decide who next.
!
!-----------------------------------------------------------------------
!
    intcnt(level) =   intcnt(level) - 1
    icheck(level)  =  icheck(level) - 1

    IF (level < lfine) level = level + 1
    GO TO 30
!
!-----------------------------------------------------------------------
!
!  One complete coarse grid integration cycle done.
!
!-----------------------------------------------------------------------
!

    110   CONTINUE

!
!-----------------------------------------------------------------------
!
!  Update all grids at levels finest-1 to level
!
!-----------------------------------------------------------------------
!
!  print*,'calling update after statement 110'
!  write(6,'(''Calling UPDATE for level '',i5)') level

    CALL update( level )

!  write(6,'(''Called  UPDATE for level '',i5)') level

    time    = time   + delt0

!
!-----------------------------------------------------------------------
!
!  Data output
!
!-----------------------------------------------------------------------
!
    IF ( ( (iout(1) > 0 .AND. MOD(nstepc,iout(1)) == 0) .OR.            &
             (iout(2) > 0 .AND. MOD(nstepc,iout(2)) == 0) .OR.          &
             (iout(3) > 0 .AND. MOD(nstepc,iout(3)) == 0) .OR.          &
             (iout(4) > 0 .AND. MOD(nstepc,iout(4)) == 0) ) .AND.       &
           nstepc >= iout(6) )THEN

      CALL usrout1( 1,lfine, 0 )

    END IF
!
!-----------------------------------------------------------------------
!
!  Run time plotting:
!
!-----------------------------------------------------------------------
!

!  IF( iout(5).gt.0 .and. mod(nstepc,iout(5)).eq.0) ) call pltall(3)

!
!-----------------------------------------------------------------------
!
!  Dump out restart data every iout(3) base grid time steps
!
!-----------------------------------------------------------------------
!

    IF ( iout(3) > 0 .AND. MOD(nstepc,iout(3)) == 0 .AND.               &
           nstepc >= iout(6) ) THEN

      WRITE(6,'(''  DUMPING RESTRART DATA AT TIME='',F10.2)') time
      CALL rstrdwr(2,time)

      WRITE(6,'(1x,a,i20)') 'Alloc dimension was     ',lstore

    END IF

!
!-----------------------------------------------------------------------
!
!  Print the cpu statistics
!
!-----------------------------------------------------------------------
!
!  cpu_simulation = f_cputime() - cpu_simulation
!  call prtcpu(1,6)

!
!-----------------------------------------------------------------------
!
!  End of one entire coarse grid time step integration cycle
!
!-----------------------------------------------------------------------
!
    WRITE(6,'(1x,a,i20)') 'Maximum alloc used was  ',ihighwater

  END DO

  RETURN
END SUBROUTINE tick
!

SUBROUTINE regrd1( lbase,rdbpts )
!
!-----------------------------------------------------------------------
!
!  Not implemented.
!
!-----------------------------------------------------------------------
!
  LOGICAL :: rdbpts
  WRITE(6,'(''  IN REGRD1 FOR LBASE '',I4)') lbase
  WRITE(6,'(''  RDBPTS IS  '',L1)') rdbpts

  RETURN
END SUBROUTINE regrd1
!

SUBROUTINE update( level )
!
!-----------------------------------------------------------------------
!
!  Update all grids at levels finest-1 to level
!
!-----------------------------------------------------------------------
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agricpu.inc'
  INCLUDE 'agricst.inc'

  cpu0 = f_cputime()

  lget = lfine - 1

  4   IF(lget < level) GO TO 999

  mptr = lstart(lget)
  5     IF (mptr == 0) GO TO 998
!
  IF(verbose6)WRITE(6,'(''  CALLING UPDGRD FOR GRID '',I5)') mptr

  CALL updgrd( mptr )

  IF(verbose6)WRITE(6,'(''  UPDGRD CALLED FOR GRID '',I5)') mptr
!
!
!-----------------------------------------------------------------------
!
!    Take care of next grid
!
!-----------------------------------------------------------------------
!
  mptr = node(10,mptr)
  GO TO 5

  998   CONTINUE
!
!-----------------------------------------------------------------------
!
!    Update next coarser level
!
!-----------------------------------------------------------------------
!
  lget = lget - 1
  GO TO 4
!
  999   CONTINUE

  cpu_update = cpu_update + f_cputime() - cpu0

  RETURN
END SUBROUTINE update
!

SUBROUTINE advanc( level )
!
!-----------------------------------------------------------------------
!
!  This routine advances all grids at level 'level' a single timestep.
!  It also calls tha appropriate boundary conditions routines for
!  interpolating fine grid boundary values
!
!  Interpolate boundary values for all grids at this level from
!  coarser grids if necessary
!
!-----------------------------------------------------------------------
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agricst.inc'
  LOGICAL :: boundc, samlvl

  boundc = .false.

  IF (level > 1) THEN

    irr = nint(possk(level-1)/possk(level))
    dtl = possk(level)
    mptr = lstart(level)
    mptrc = lstart(level-1)
    istp = nint((rnode(20,mptrc)-rnode(20,mptr))/dtl)

    tmwght = 1./(1.+intratt)  ! Time interpolation coefficient
!
!-----------------------------------------------------------------------
!
!  Set boundc to .true. for the first of a number of fine grid time
!  steps that will make one coarse grid time step.
!
!  Only for this first time step, boundary conditions are interpolated
!  from the coarse grid.
!
!-----------------------------------------------------------------------
!
    IF(istp == irr) boundc = .true.

    IF( boundc .AND. verbose6 ) THEN
      WRITE(6,*)' tmwght=',tmwght,                                      &
                ' level=', level,' intratt= ',intratt
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Loop through all grids at level 'level'
!
!-----------------------------------------------------------------------
!
  mptr = lstart(level)

  15    CONTINUE

  IF( mptr /= 0) THEN
!
!-----------------------------------------------------------------------
!
!  Call coarse grid boundary interpolation routine to obtain the
!  boundary time tendencies for the fine grid when boundc=.true.
!
!-----------------------------------------------------------------------
!
    IF( boundc ) THEN

      IF( verbose6 ) WRITE(6,'('' CALLING BOUNDARY INTERP ROUTINE FOR '', &
      &           '' interp from coarser grids '')')

      samlvl = .false.   ! For grids at different levels.

      CALL updbc( mptr, samlvl, tmwght )

    END IF
!
!-----------------------------------------------------------------------
!
!  Call solver for mptr
!
!-----------------------------------------------------------------------
!
    IF(verbose6)WRITE(6,'('' CALLING ARPSOLVE FOR GRID '',I4)')mptr

    CALL arpsolve( mptr )
!
!-----------------------------------------------------------------------
!
!  Do the next grid at this level
!
!-----------------------------------------------------------------------
!
    mptr = node(10,mptr)
    GO TO 15

  END IF
!
!-----------------------------------------------------------------------
!
!  Exchange boundary values between all grids at this level
!  if necessary
!
!-----------------------------------------------------------------------
!
  CALL exchng( level )
!
  RETURN
END SUBROUTINE advanc
