SUBROUTINE fillngrd(levl,levd)
!
!  this routine interpolates initial fields for the
!  new grids at level levl
!

  INTEGER :: levl,levd

  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agricst.inc'

  LOGICAL :: sourced
!
!  interp for all grids
!
  IF( verbose6 ) WRITE(6,'('' IN INITGR, NEWSTL(LEVL) = '',I4)')newstl(levl)
  WRITE(6,'('' LEVL = '',I4)') levl
!
  mptr = newstl(levl)
  5    IF (mptr == 0) GO TO 999
!
!  set storage for the new grid
!
  CALL resett
  CALL setstr( mptr )
  CALL settmp
!
!  first interp off of the coarse grids,
!  then from grids at same level
!
!cc     do 300 levint=max0(1,levl-2),levl
!
!?????
!
  DO levint=MAX0(1,levl-1),levl


    sourced = .false.
    IF( levint >= levd ) sourced = .true.
    msrc = lstart(levint)
    10       IF( msrc == 0 ) CYCLE
!
    IF( verbose6 ) WRITE(6,*) ' new fine grid interp ',mptr,msrc
!
!  interp from grid msrc
!
    IF ( sourced ) CALL addgrd( msrc )
    CALL filgrd( mptr,msrc )
    IF ( sourced ) CALL dmpgrd( msrc,.false. )
!
!  check for next source grid on level
!
    msrc = node(10,msrc)
    GO TO 10
!
  END DO
!
!  take care of next grid
!
  CALL dmpgrd( mptr,.true. )
  mptr = node(10,mptr)
  GO TO 5
  999    CONTINUE
!
!  now set all the constant arrays for the new grid
!
  mptr = newstl(levl)
  105     IF (mptr == 0) GO TO 9999
  CALL addgrd( mptr )
  CALL flcnst( mptr )
  mptr = node(10,mptr)
  GO TO 105
  9999    CONTINUE
!
! we're finished here
!
! yuhe: level levl has been filled. Set lstart(levl) in order to use
!    this level to source finer levels.
!
  lstart(levl) = newstl(levl)

  RETURN
END SUBROUTINE fillngrd
