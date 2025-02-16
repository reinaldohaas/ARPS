!
! --------------------------------------------------------------------
!

  FUNCTION igetsp (nwords)
!
  INCLUDE 'manage.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agricst.inc'
!
!  allocate contiguous space of length nword in main storage array
!  alloc. user code (or pointer to the owner of this storage)
!  is  mptr.  lenf = current length of lfree list.
!
!
!  find first fit from free space list
!
  itake = 0
  DO i = 1, lenf
    IF (lfree(i,2) < nwords) CYCLE
    itake = i
    GO TO 25
  END DO
  GO TO 900
!
!  put remaining words back on list.
!
  25   left = lfree(itake,2) - nwords
  igetsp = lfree(itake,1)
!
  IF (left <= 0) GO TO 30
  lfree(itake,2) = left
  lfree(itake,1) = lfree(itake,1) + nwords
  GO TO 99
!
!  item is totally removed.  move next items in list up one.
!
  30   lenf = lenf - 1
  DO i = itake, lenf
    lfree(i,1) = lfree(i+1,1)
    lfree(i,2) = lfree(i+1,2)
  END DO
  GO TO 99
!
  900  WRITE(6,901) nwords
  901  FORMAT('  require ',i16,' words - either none left or not big',  &
                '  enough space')
  WRITE(6,902) ((lfree(i,j),j=1,2),i=1,lenf)
  902  FORMAT(' free list: ',//,2X,50(i10,4X,i10,/,2X))
  WRITE(6,*) ' ******* you need more space !!!! ********** '
  WRITE(6,*) ' increase lstore in agrialloc.inc '
!c      call outtre(mstart,.false.,nvar)
  STOP
!
  99   IF (verbose4) WRITE(6,100) nwords, igetsp
  ihighwater = MAX0( ihighwater,lfree(itake,1) )
  100  FORMAT('  allocating ',i9,' words in location ',i9)
!   if (verbose4) then
!      print*, 'ending call to igetsp'
!      print*, 'nwords in call = ', nwords
!      print*, 'lenf at present = ', lenf
!      print *, 'ihighwater =  ', ihighwater
!      do 666 i=1, lenf
!         print *, 'i, lfree(i,1), lfree(i,2)',i,lfree(i,1),lfree(i,2)
! 666     end do
!   end if
  RETURN
  END FUNCTION igetsp
!
! --------------------------------------------------------------------
!

SUBROUTINE reclam (INDEX, nwords)
!
!  return of space. add to free list.
!  iplace points to next item on free list with larger index than
!  the item reclaiming, unless said item is greater then
!  everything on the list.
!
  INCLUDE 'manage.inc'
!
  10   DO i = 1, lenf
    iplace  = i
    IF (lfree(i,1) > INDEX) GO TO 30
  END DO
  WRITE(6,902)
  902     FORMAT(' no insertion pointer into freelist. error stop')
  PRINT *, 'failed iptr= ', INDEX
  PRINT *, 'lenf =', lenf
  PRINT *, 'lfree(lenf,1) = ', lfree(lenf,1)
  PRINT *, 'nwords  = ', nwords
  STOP
!
!  check previous segment for merging
!
  30      iprev = iplace - 1
  IF (lfree(iprev,1)+lfree(iprev,2) < INDEX) GO TO 40
  lfree(iprev,2) = lfree(iprev,2) + nwords
  GO TO 50
!
!  check after segment - no previous merge case
!
  40   nexti = INDEX + nwords
  IF (lfree(iplace,1) /= nexti) GO TO 70
  lfree(iplace,1) = INDEX
  lfree(iplace,2) = lfree(iplace,2) + nwords
  GO TO 99
!
!  check following segment - yes previous merge case
!
  50   nexti = INDEX + nwords
  IF (lfree(iplace,1) /= nexti) GO TO 99
!
! forward merge as well, bump all down 1
!
  lfree(iprev,2) = lfree(iprev,2)+lfree(iplace,2)
  ipp1           = iplace + 1
  DO i = ipp1, lenf
    lfree(i-1,1) = lfree(i,1)
    lfree(i-1,2) = lfree(i,2)
  END DO
  lenf = lenf - 1
  GO TO 99
!
!  no merges case - insert and bump future segments up to make room
!

  70   IF (lenf == idimf) GO TO 900
  DO ii = iplace, lenf
    i          = lenf + 1 - ii + iplace
    lfree(i,1) = lfree(i-1,1)
    lfree(i,2) = lfree(i-1,2)
  END DO
  lenf            = lenf + 1
  lfree(iplace,1) = INDEX
  lfree(iplace,2) = nwords
  GO TO 99
!
  900  WRITE(6,901) idimf
  901  FORMAT('  free list full with ',i5,' items')
  STOP
!
  99   IF (sprint) WRITE(6,100) nwords, INDEX
  100  FORMAT('     reclaiming ',i9,' words at loc. ',i9)
  RETURN
END SUBROUTINE reclam
!
! --------------------------------------------------------------------
!

SUBROUTINE resett
  INCLUDE 'manage.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agricst.inc'
!
!  reset the temp space
!
  IF( verbose4 ) WRITE(6,'(''  RESETTING TEMP SPACE, NTEMP= '',I8)')ntemp
  DO  i  = 1, idimf
    lfree(i,1) = 0
    lfree(i,2) = 0
  END DO
!
  lfree(3,1) = + 1
  lfree(2,1) = ntemp
  lfree(2,2) = lstore - ntemp + 1
  lenf       = 3
!
  RETURN
END SUBROUTINE resett
!
! --------------------------------------------------------------------
!

SUBROUTINE stst1
!
  INCLUDE 'manage.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
!
!  intialize a few variables needed before calling user set up
!  routine domain. finish in stst2 after.
!  the spatial and temporal stepsizes are set. the node array
!  is kept as a linked list of free nodes.  "ndfree" points to the
!  head of the list, i.e.-first free node.  use first row of each
!  col to hold this pointer, set by the macro "nextfree".
!  the free space list, managed in lfree, will have first and
!  last positions filled with an allocation of zero words,
!  to avoid boundary cases.
!
  ndfree = 1
  DO i   = 1,30
    node(2,i) = i+1
  END DO
!
! the last free node will have a null pointer

  node(2,30) = 0
!
!  initialize linked list of alloc storage as well.
!  first and last locations are dummy placeholders of zero words
!  of allocation each, to avoid boundary cases.
!
  idimf = 50
  DO  i  = 1, idimf
    lfree(i,1) = 0
    lfree(i,2) = 0
  END DO
!
  lfree(3,1) = + 1
  lfree(2,1) = 1
  lfree(2,2) = lstore
  lenf       = 3
!
  RETURN
END SUBROUTINE stst1
!
! --------------------------------------------------------------------
!

SUBROUTINE stst2
!
  INCLUDE 'manage.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
!
! stst = finish initializing spatial and counting arrays.
!
  levmlt(1)  =  1
  lev        = 2
  rr         = FLOAT(intrat)
  rrt        = FLOAT(intratt)
  10   IF (lev > mxnest) GO TO 20
  hxposs(lev) = hxposs(lev-1) / rr
  hyposs(lev) = hyposs(lev-1) / rr
  possk (lev) = possk (lev-1) / rrt
  levmlt(lev) = intratt
  lstart(lev) = 0
  lback(lev)  = 0
  newstl(lev) = 0
  lev         = lev + 1
  GO TO 10
  20   CONTINUE
!
! after kcheck integrations of parent grid, move its refinements.
! finest level grid never needs to have its finer subgrids moved.
!
  DO i = 1,3
    ncheck(i) = kcheck
  END DO
  ncheck(mxnest) = 1000000000
!
  RETURN
END SUBROUTINE stst2
!
! --------------------------------------------------------------------
!

SUBROUTINE stst3
!
  INCLUDE 'manage.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
!
! stst = finish initializing spatial and counting arrays.
!
! this routine is also called with each restart.  In this way
! the temporal refinement can be changed on the fly.
! The spatial refinement cannot be changed on the fly!!!!!
!
  levmlt(1)  =  1
  lev        = 2
  rr         = FLOAT(intrat)
  rrt        = FLOAT(intratt)
  10   IF (lev > mxnest) GO TO 20
  hxposs(lev) = hxposs(lev-1) / rr
  hyposs(lev) = hyposs(lev-1) / rr
  possk (lev) = possk (lev-1) / rrt
  levmlt(lev) = intratt
  lev         = lev + 1
  GO TO 10
  20   CONTINUE
!
! after kcheck integrations of parent grid, move its refinements.
! finest level grid never needs to have its finer subgrids moved.
!
  DO i = 1,3
    ncheck(i) = kcheck
  END DO
  ncheck(mxnest) = 1000000000
!
  RETURN
END SUBROUTINE stst3
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION nodget(dummy)
!
  INCLUDE 'nodal.inc'
!
! nodget =  get first free node of the linked list kept in node
!         array. adjust pointers accordingly.
!
  IF (ndfree /= 0) GO TO 10
  WRITE(6,100)
  100       FORMAT(20H  out of nodal SPACE)
  STOP
!
  10     nodget         = ndfree
  ndfree         = node(2,ndfree)
!
!  initialize nodal block
!
  DO i        = 1,25
    node(i,nodget) = 0
  END DO
  DO i         = 1,35
    rnode(i,nodget) = 0.0
  END DO
!
  RETURN
  END FUNCTION nodget
!
! --------------------------------------------------------------------
!

SUBROUTINE putnod (mptr)
!
  INCLUDE 'nodal.inc'
!
! putnod = return mptr node to the linked list kept in node array.
!
  node(2, mptr) = ndfree
  ndfree        = mptr
!
  RETURN
END SUBROUTINE putnod
!
