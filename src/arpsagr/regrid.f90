SUBROUTINE regrid(lbase,rdbpts)
!
  LOGICAL :: rdbpts
  INCLUDE 'nodal.inc'
  INCLUDE 'agricpu.inc'
  INCLUDE 'agricst.inc'
!
! ******************************************************************
!  REGRID = FLAG POINTS ON EACH GRID WITH A LEVEL > = LBASE.
!  CLUSTER THEM, AND FIT NEW SUBGRIDS AROUND THE CLUSTERS.
!  THE LBASE GRIDS STAY FIXED DURING REGRIDDING OPERATION.
!  WHEN A PARENT GRID HAS ITS ERROR ESTIMATED, ADD ITS KID GRID
!  INFORMATION TO THE ERROR GRID BEFORE CLUSTERING. (PROJECT)
!  ORDER OF GRID EXAMINATION - ALL GRIDS AT THE SAME LEVEL, THEN
!  DO THE NEXT COARSER LEVEL.
!
!  NOTE:  THE ACTUAL GRIDFITTING IS DONE BY GRDFIT.  IN THIS
!      VERSION WE DO NOT HAVE POINT FLAGGING CAPABILITY IN
!      3-D (NOR GRIDFITTING CAPABILITY IN 3-D) HENCE THE
!      NEW GRIDS ARE SPECIFIED BY THE USER.  SEE GRDFIT.
!
!      RDBPTS = true for reading badpoints from file
!
! LOCAL VARIABLES:
!  LCHECK = THE LEVEL BEING EXAMINED.
!  LFNEW  = FINEST GRID TO BE. WILL REPLACE LFINE.
! GLOBAL
!    MSTART  = START OF VERY COARSEST GRIDS.
! *****************************************************************
!
  cpu0 = f_cputime()
!
!
  PRINT*,'REGRID called with LBASE=', lbase,' mxnest=',mxnest
  PRINT*,'LFINE =', lfine

!
! yuhe: LCHECK should go from lower to higher so that the grid numbers
!    will be set from lower to higher.
!
!  LCHECK    = MIN0(LFINE,MXNEST-1)
  lcheck    = lbase
  lfnew     = lbase

  DO i   = 1, mxnest
    newstl(i) = 0
  END DO

  time      = rnode(20, lstart(lbase))

  20   CONTINUE

  PRINT*,' lcheck =',lcheck,' lbase= ',lbase

!  IF (LCHECK .ge. LBASE) THEN
  IF ( lcheck <= MIN0(lfine,mxnest-1) ) THEN
!
! yuhe: Since lstart(lcheck) has not been set, the following statement
!    may run into a trouble.
!
!     IF (ABS(TIME-RNODE(20, LSTART(LCHECK))).gt..0001) THEN
!
!       WRITE(6,100) TIME
!100       FORMAT(41H  GRIDS AT DIFFERENT TIMES DURING REGRID ,E12.7)
!       STOP
!
!     ENDIF

    WRITE(6,'(''  IN REGRID, CALLING GRDFIT '')')

    CALL grdfit(lbase,lcheck,lfnew,time,rdbpts)

    PRINT*,' LBASE,LCHECK,LFNEW,TIME,RDBPTS =',                         &
             lbase,lcheck,lfnew,time,rdbpts, newstl(lcheck+1)

    IF (newstl(lcheck+1) /= 0) lfnew = MAX0(lcheck + 1,lfnew)
!
! note: this line is added for testing purpose
!
!    LFNEW = MAX0(LCHECK + 1,LFNEW)

!    LCHECK = LCHECK - 1
    lcheck = lcheck + 1

    GO TO 20

  END IF

  PRINT*,' LBASE,LCHECK,LFNEW,TIME,RDBPTS =',                           &
           lbase,lcheck,lfnew,time,rdbpts

  PRINT*,' END OF LEVEL LOOP '

!
!  END OF LEVEL LOOP
!
!  REMAINING TASKS LEFT IN REGRIDDING:
!
!  INTERPOLATE STORAGE FOR THE NEW GRIDS.  THE STARTING POINTERS
!  FOR EACH LEVEL ARE IN NEWSTL.
!
!  first we need to push the old data out temporarily
!
  PRINT*,' lfnew ', lfnew

  levd = lfnew
  DO levn=lfnew,lbase+1,-1

    PRINT*,' levn =', levn

    IF(lstart(levn) /= 0) THEN

      PRINT*,'calling dmplvl '

      CALL dmplvl(levn)
      levd = levn
    END IF
  END DO
!
!  now fill the new grids
!
  DO levn=lbase+1,lfnew
    IF( verbose6 ) WRITE(6,'(''  CALLING FILL FOR LEVEL '',I5)') levn

    PRINT*,' levd =', levd

    CALL fillngrd( levn, levd )
  END DO
!
!  clean get rid of disk files for the temp grid dump
!  performed by dmplvl
!
  PRINT*,' lfnew, lbase+1 ', lfnew, lbase+1

  DO levn=lfnew,lbase+1,-1

    PRINT*,' levn =',levn, lfnew, lbase+1

    mptr = lstart(levn)
    CALL cleanf( mptr )
    mptr = newstl(levn)
    CALL cleanf( mptr )
  END DO
!
!  MERGE DATA STRUCTURES (NEWSTL AND LSTART, LBACK)
!  FINISH STORAGE ALLOCATION, RECLAIM SPACE,
!  SET UP NEW SOURCE ARRAYS FOR INTERP AND UPDATE
!
  PRINT*,'calling join with LBASE, LFNEW=',lbase, lfnew

  CALL join(lbase, lfnew)
!
  cpu_regrid = cpu_regrid + f_cputime() - cpu0
!
  RETURN
END SUBROUTINE regrid
!

SUBROUTINE grdfit (lbase,lcheck,lfnew,time,rdbpts)
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agricst.inc'
  LOGICAL :: rdbpts
  INTEGER :: lcheck,prvptr
  REAL :: ugrids(11,10)
!
!  GRDFIT CALLED BY SETGRD AND REGRID TO ACTUALLY FIT THE NEW GRIDS
!      ON EACH LEVEL. LCHECK IS THE LEVEL BEING ERROR ESTIMATED
!      SO THAT LCHECK+1 WILL BE THE LEVEL OF THE NEW GRIDS.
!
!  IN THIS VERSION WE JUST CALL A USER WRITTEN ROUTINE WHICH RETURNS
!  THE NEW GRID LOCATION IN UGRID.  UGRID(1-8,*) ARE THE FOUR CORNERS
!  OF THE NEW GRID, 9 AND 10 ARE THE BOTTOM AND TOP, 11 IS DZ.
!
  lckp1 = lcheck+1

  IF(ngrdnew(lcheck) == 0) GO TO 99

  CALL rdgrds(ugrids,lckp1)
!
!  SET UP GRID STRUCTURE FOR EACH NEW GRID
!
  levnew = lckp1
  prvptr=0
  DO ig=1,ngrdnew(lckp1-1)
    mnew = nodget(dummy)
    IF(prvptr /= 0) node(10,prvptr) = mnew
    IF(prvptr == 0) newstl(levnew)  = mnew
    node(12,mnew) = prvptr
    DO ip=1,8
      rnode(ip,mnew) = ugrids(ip,ig)
    END DO
!    CALL MOMENT( RNODE(1,MNEW),UGRIDS(1,IG),4,
!    *               USAGE,0.,0.,HXPOSS(LEVNEW),HYPOSS(LEVNEW) )
    rnode(20,mnew) = time
    rnode(29,mnew) = ugrids(9,ig)
    rnode(30,mnew) = ugrids(10,ig)
    rnode(31,mnew) = ugrids(11,ig)
    node(4,mnew) = levnew
    node(9,mnew) = 1
    prvptr = mnew
  END DO
!
  CALL sethk (newstl(levnew))
  CALL setrot(newstl(levnew))
!
!
  99   RETURN
END SUBROUTINE grdfit
!

SUBROUTINE rdgrds(ug,ln)
  DIMENSION ug(11,10), tmpg(8)
  INCLUDE 'nodal.inc'
  INCLUDE 'agricst.inc'
!
!  this routine returns the corners for the new grids
!  ln is the level for the new grids.  ngrdnew is the
!  numer of new grids on ln passed through common block
!  agri100.  ug(1-8) contains the x-y coordinates of the
!  corners, 9-10 the bottom and top height and 11 is dz.
!
!
!  we will read in the new grid information at present
!  xc and yc is the center point, xl and yl the
!  axis lengths, angle is the rotation angle, zo and
!  zt are the bottom and top of the grid boundaries
!
!  these are in coarse-grid normalized coordinates.
!
!
  WRITE(6,*) '  regridding'
  DO ig=1,ngrdnew(ln-1)

    pio2 = 1.570796
    angle = -pio2*gangle(ig,ln-1)/90.
    cosang = COS(angle)
    sinang = -SIN(angle)
    xld2 = 0.5*ixln(ig,ln-1)
    yld2 = 0.5*jyln(ig,ln-1)

    ug(1,ig) = (ixc(ig,ln-1) - xld2*cosang + yld2*sinang -1 ) *         &
                                   rnode( 9,mstart)
    ug(2,ig) = (jyc(ig,ln-1) - xld2*sinang - yld2*cosang -1 ) *         &
                                   rnode( 9,mstart)

    ug(3,ig) = (ixc(ig,ln-1) - xld2*cosang - yld2*sinang -1 ) *         &
                                   rnode( 9,mstart)
    ug(4,ig) = (jyc(ig,ln-1) - xld2*sinang + yld2*cosang -1 ) *         &
                                   rnode( 9,mstart)

    ug(5,ig) = (ixc(ig,ln-1) + xld2*cosang - yld2*sinang -1 ) *         &
                                   rnode( 9,mstart)
    ug(6,ig) = (jyc(ig,ln-1) + xld2*sinang + yld2*cosang -1 ) *         &
                                   rnode( 9,mstart)

    ug(7,ig) = (ixc(ig,ln-1) + xld2*cosang + yld2*sinang -1 ) *         &
                                   rnode( 9,mstart)
    ug(8,ig) = (jyc(ig,ln-1) + xld2*sinang - yld2*cosang -1 ) *         &
                                   rnode( 9,mstart)
!
! here we redo the points such that the point with max x is
! always the first point.  this means that the grids will
! always have an angle between 0 and -90 degrees in the interface.
! this does not change the actual grid placement in any way.
!
    i1 = 1
    IF(ug(3,ig) < ug(i1,ig)) i1 = 3
    IF(ug(5,ig) < ug(i1,ig)) i1 = 5
    IF(ug(7,ig) < ug(i1,ig)) i1 = 7
!
    IF(i1 > 1) THEN
      DO  i=i1,i1+7
        in = i
        IF(i > 8) in = i-8
        ii = i-i1+1
        tmpg(ii) = ug(in,ig)
      END DO
      DO  i=1,8
        ug(i,ig) = tmpg(i)
      END DO
    END IF
!
    ug(9,ig) = 0.
    ug(10,ig)= rnode(30,mstart)
    ug(11,ig)= rnode(31,mstart)
    WRITE(6,11)ig,(ug(ii,ig),ii=1,11)
  END DO

11  FORMAT(2X,' new grid ',i4,', location information ',3(/,2X,4E12.5) )
  RETURN
END SUBROUTINE rdgrds
!

SUBROUTINE join(lbase, lfnew)
!
  INCLUDE 'nodal.inc'
  LOGICAL :: printl
  INTEGER :: ptr
!
!  JOIN = HOOK UP OLD AND NEW LEVEL DATA STRUCTURES.
!      SET PREVLEVEL POINTERS.
!      RESET GLOBAL VARIABLE  LFINE.
!      RETURN FILES FROM OLD GRIDS.
!
!
  l = lbase + 1
  10   IF (l > lfine) GO TO 40
  mptr = lstart(l)
  20       IF (mptr == 0) GO TO 30
  mold = mptr
  mptr = node(10, mptr)
!          CALL PUTNOD(MOLD)
  GO TO 20
  30     CONTINUE
  l = l + 1
  GO TO 10
  40   CONTINUE
!
! MERGE NEWSTL INTO LSTART FOR FINISHING TOUCHES OF GRID STRUCTURE
!
  l = lbase + 1
  50   IF (l > 10) GO TO 60
  lstart(l) = newstl(l)
  l = l + 1
  GO TO 50
  60   CONTINUE

  lfine = lfnew
!
!  GRID STRUCTURE NOW COMPLETE AGAIN. SAFE TO PRINT, ETC. ASSUMING
!  THINGS INITIALIZED TO ZERO IN NODGET.
!
! here we'd set whatever else needs to be set
!
  printl = .true.
  l = lbase + 1
  70   IF (l > lfine) GO TO 105
  mptr = newstl(l)
  80       IF (mptr == 0) GO TO 90
  lm1=l-1
!           CALL SRCLST( MPTR,MPTRSC,DXYSRC,NUMBS,LM1,L )
!           CALL ZSOURCE( IZSRC(1,1,MPTR),MPTRSC,MPTR,NUMBS,
!  *                      BOUNDI,PRINTL )
  mptr = node(10, mptr)
  GO TO 80
  90       CONTINUE
  l = l + 1
  GO TO 70
  105   CONTINUE
!
! SET THE LBACK POINTER
!
  100  l      = lbase+1
  110  IF (l > lfine) GO TO 130
  ptr    = lstart(l)
  120      CONTINUE
  lback(l) = ptr
  ptr      = node(10,ptr)
  IF (ptr /=  0) GO TO 120
  l     = l + 1
  GO TO 110
  130  CONTINUE
!
!
  RETURN
END SUBROUTINE join
!

SUBROUTINE sethk(msave)
  INCLUDE 'nodal.inc'
  INTEGER :: mptr
  INTEGER :: pptr,n,mlevel
  REAL :: inner,norm1,norm2,cossqu,h1,h2
!  LOGICAL  OVRLAP
! ******************************************************************
! SETHK  - SET THE H'S, K, AND MAXNUMROW, COL FOR THE GRIDS
!       STARTING AT MSAVE.
! INPUT PARAMETER:
!    MPTR   - STARTING PTR TO GRID TO FIND HROW, HCOL, K, MAXI, MAXJ FOR
!
! NOTE ON ALGORITHM:
! SINCE THE STEP SIZES MAY BE DIFFERENT DIRECTIONS, A SMOOTH CHANGE
! FROM 0 TO PI/2 ORIENTATION ANGLE IS MADE BY COMPUTING THE ANGLE
! BETWEEN MPTR AND AN (OLD) PARENT.  CALL THIS ANGLE ALPHA.
! WE THEN FIND WHICH OF THE H'S BY INSISTING THAT IF ALPHA = 0,
! THEN THE NEW H IS JUST H(X OR Y)POSS; IF ALPHA = PI/2, THEN THE
! NEW H IS JUST H(Y OR X)POSS; AND THAT THE CHANGE IS SMOOTH.
! THE CHOICE MADE WAS (NOTE THAT THIS HAS NOTHING TO DUE WITH
! ROTATIONS):
!    NEW H1 = H1*COS**2(ALPHA) + H2*SIN**2(ALPHA)
! AND SIMILIARLY FOR H2.
! HERE THE H1 AND H2 ARE THE "APPROPRIATE" CHOICE OF H(X OR Y)POSS
! (BASED ON THE ORIENTATION OF MPTR)
! SPECIAL NOTE:
!    IT IS USUALLY THE CASE THAT THE SIDES OF THE RECTANGLE ARE NOT
! AN INTEGRAL MULTIPLE OF THE STEP SIZES.  SETHK INCREASES THE H'S
! TO THE SMALLEST H > ORIGINAL H. THIS IS TO KEEP THE METHOD STABLE,
! WHICH WOULD NOT BE TRUE IF WE USED AN H < ORIGINAL H.
! ******************************************************************
  mlevel = node(4,msave)
  mptr   = msave
!
  23000 IF(.NOT. (mptr /= 0))GO TO 23001
  pptr   = lstart(1)
!23002     IF(.NOT. (.NOT. OVRLAP(PPTR,MPTR)))GO TO 23003
!       PPTR = NODE(10, PPTR)
!       GO TO 23002
!23003     CONTINUE
!
!       COMPUTE INNER PRODUCT AND LENGTHS OF SIDES
  inner  = ( rnode(7,pptr)-rnode(1,pptr) ) *( rnode(7,mptr)-            &
      rnode(1,mptr) ) +( rnode(8,pptr)-rnode(2,pptr) ) *( rnode(8,mptr) &
      -rnode(2,mptr) )
  norm1  = ( ( rnode(7,pptr)-rnode(1,pptr) ) ** 2 +( rnode(8,           &
      pptr)-rnode(2,pptr) ) ** 2 )
  norm2  = ( ( rnode(7,mptr)-rnode(1,mptr) ) ** 2 +( rnode(8,           &
      mptr)-rnode(2,mptr) ) ** 2 )
!
  cossqu = (inner/norm1) * (inner/norm2)
!
  IF (cossqu < 0.5) GO TO 23004
  h1     = hxposs(mlevel)
  h2     = hyposs(mlevel)
  GO TO 23005
  23004     CONTINUE
  h1     = hyposs(mlevel)
  h2     = hxposs(mlevel)
  23005     CONTINUE
!
  hx    = (h1-h2)*cossqu + h2
  hy    = (h2-h1)*cossqu + h1
  rnode(11,mptr) = possk(mlevel)
!
!       NOW, WE ADJUST THE H'S SO THAT N*H = LENGTH OF SIDE FOR SOME
!       INTEGER N
!
  norm2 = SQRT( (rnode(3,mptr)-rnode(1,mptr))**2 +(rnode(4,mptr         &
      )-rnode(2,mptr))**2 )
  norm1 = SQRT( (rnode(7,mptr)-rnode(1,mptr))**2 +(rnode(8,mptr         &
      )-rnode(2,mptr))**2 )
!
  n     = IFIX(norm1/hx+.01 + 1)
!
! SET N,M ODD FOR ERREST (2DX GRID OVERLY)
!
!   IF(MOD(N,2) .NE. 0) N=N-1
  rnode(9,mptr) = norm1 / FLOAT(n-1)
  node(5,mptr) = n
!
  n     = IFIX(norm2/hy+.01 + 1)
!   IF(MOD(N,2) .NE. 0) N=N-1
  rnode(10,mptr) = norm2 / FLOAT(n-1)
  node(6,mptr) = n
!
!  NOW SET NUMBER OF LAYERS FOR THIS GRID
!
  node(14,mptr) = nint( (rnode(30,mptr)-rnode(29,mptr))                 &
                          /rnode(31,mptr) ) + 1
!
  mptr = node(10,mptr)
  GO TO 23000
  23001 CONTINUE
  RETURN
END SUBROUTINE sethk
!
