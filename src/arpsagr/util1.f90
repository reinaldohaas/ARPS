SUBROUTINE moment (rect,badpts,npt,usage,rotate,buffer,hx,hy)
!
  REAL :: rect(28),badpts(2,npt),evec1(2),evec2(2),buffer,              &
           rotate
!
! input parameters:
!  badpts      = x,y coords of flagged badpts grouped into clusters
!                are in the first two rows
!  npt         = num. of badpts. in the cluster.
!  buffer      = expand new grids with buffer zone of size buffer
!                this prevents degenerate point rectangles.
! output parameters:
!  usage       = ratio of flagged to unflagged badpts. in new grid
!                measures goodness of fit and clustering
!    rect( )      = stores some info. for grid created herein.
!                sometimes rect = rnode, sometimes = temp. array.
!                depending on calling prog. (grdfit or expand)
!
!
!  moment = compute 2nd moments about the mean of flagged badpts.
! fit an ellipse. the eigenvectors determine rectangle orientation.
! same some info., even tho. usage might be low and rect. scrapped.
!
  rn = FLOAT(npt)
  xmean = 0.0
  ymean = 0.0
  cross = 0.0
  xsq   = 0.0
  ysq   = 0.0
  DO ipt = 1, npt
    xmean = xmean + badpts(1,ipt)
    ymean = ymean + badpts(2,ipt)
    xsq   = xsq   + badpts(1,ipt)**2
    ysq   = ysq   + badpts(2,ipt)**2
    cross = cross + badpts(1,ipt)*badpts(2,ipt)
  END DO
!
  xmean = xmean / rn
  ymean = ymean / rn
  xm2   = xmean ** 2
  ym2   = ymean ** 2
  xsq   = xsq / rn
  ysq   = ysq / rn
  cross = cross / rn - xmean * ymean
!
  tl    = xsq - xm2
  br    = ysq - ym2
!
  eval1 = ( (tl+br) + SQRT((tl-br)**2 + 4.0*cross**2) ) / 2.0
  eval2 = ( (tl+br) - SQRT((tl-br)**2 + 4.0*cross**2) ) / 2.0
!
! matrix = <tl, cross;  cross, br>.
! determine eigenvectors.  normalize for ease in later computation.
! make 1st evector lie in southeast quadrant.
! make 2nd evector lie in northeast quadrant.
!
  IF (ABS(cross) <= rotate) GO TO 30
  evec1(2) = (eval1-tl) / cross
  evec2(2) = (eval2-tl) / cross
  IF (evec1(2) <= evec2(2)) GO TO 20
  temp     = evec1(2)
  evec1(2) = evec2(2)
  evec2(2) = temp
  20       CONTINUE
  fac1     = SQRT(1.0+evec1(2)**2)
  evec1(1) = 1.0 / fac1
  evec1(2) = evec1(2) / fac1
  fac2     = SQRT(1.0 + evec2(2)**2)
  evec2(1) = 1.0 / fac2
  evec2(2) = evec2(2) / fac2
  GO TO 40
  30   CONTINUE
  evec1(1) = 1.0
  evec2(1) = 0.0
  evec1(2) = 0.0
  evec2(2) = 1.0
  40   CONTINUE
!
! compute length of enclosing rectangles to include all flagged badpts.
! take projection (dotproduct) of bad pt. on evectors to do this.
! expand on all sides by 'buffer' zone.  this also prevents
!  case of degenerate rectangle (line).
!
  emx1 = badpts(1,1)*evec1(1) + badpts(2,1)*evec1(2)
  emn1 = emx1
  emx2 = badpts(1,1)*evec2(1) + badpts(2,1)*evec2(2)
  emn2 = emx2
  IF (npt == 1) GO TO 90
  DO ipt = 2, npt
    dot1  = badpts(1,ipt)*evec1(1) + badpts(2,ipt)*evec1(2)
    dot2  = badpts(1,ipt)*evec2(1) + badpts(2,ipt)*evec2(2)
    IF (dot1 <= emx1) GO TO 50
    emx1 = dot1
    50       IF (dot1 >= emn1) GO TO 60
    emn1 = dot1
    60       IF (dot2 <= emx2) GO TO 70
    emx2 = dot2
    70       IF (dot2 >= emn2) CYCLE
    emn2 = dot2
  END DO
  90   emx1  = emx1 + buffer
  emx2  = emx2 + buffer
  emn1  = emn1 - buffer
  emn2  = emn2 - buffer
!
! from length of the sides, determine the 4 rect. corners
! number the corners clockwise
!
  rect(1) = emn1*evec1(1)+emn2*evec2(1)
  rect(2) = emn1*evec1(2)+emn2*evec2(2)
  rect(3) = emn1*evec1(1)+emx2*evec2(1)
  rect(4) = emn1*evec1(2)+emx2*evec2(2)
  rect(5) = emx1*evec1(1)+emx2*evec2(1)
  rect(6) = emx1*evec1(2)+emx2*evec2(2)
  rect(7) = emx1*evec1(1)+emn2*evec2(1)
  rect(8) = emx1*evec1(2)+emn2*evec2(2)
!
!  set some other fields here that would be difficult elsewhere.
!
  rect(16)  = emx1  +.00001
  rect(17)  = emn1  -.00001
  rect(18)  = emx2  +.00001
  rect(19)  = emn2  -.00001
  rect(12)  = evec1(1)
  rect(13)  = evec1(2)
  rect(14) = evec2(1)
  rect(15) = evec2(2)
  dist = SQRT( (rect(7)-rect(1))**2                                     &
              +(rect(8)-rect(2))**2 )
  rect( 21)=(rect(7)-rect(1))/dist
  rect( 22)=(rect(8)-rect(2))/dist
!
! compute volume cutoff ratio
!
  iside1 = (hx+emx1-emn1-2.0*buffer+.00001) / hx
  iside2 = (hy+emx2-emn2-2.0*buffer+.00001) / hy
  gpall  = iside1 * iside2
  usage  = rn / gpall
!
  RETURN
END SUBROUTINE moment
!
! --------------------------------------------------------------------
!

SUBROUTINE outtre(mlev,ddsc,dmpc)
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agricpu.inc'
  INCLUDE 'agricst.inc'
  LOGICAL :: ddsc,dmpc
!
! ******************************************************************
! OUTTRE - OUTPUT SUBTREE
! INPUT PARAMETERS:
!    MLEV   - PTR TO SUBTREE TO OUTPUT I.E., START AT LEVEL(MLEV)
! TREE IS OUTPUT FROM 'LEVEL' TO FINEST LEVEL.
! ******************************************************************
!
!
!  dump out generic grid data if asked for
!
  IF(ddsc) THEN
    CALL dmpdsc
    CALL dmpstr
  END IF
!
! dump out tree if needed
!
  IF( .NOT. verbose1 ) RETURN
!
  time = rnode(20,mlev)/60.

  WRITE(6,1) time
  1     FORMAT(//,1X,'  GRID STRUCTURE AT  ',f10.3,' MINUTES ')
!
  level = node(4, mlev)
!
  IF( verbose6 ) WRITE(6,'(''  LEVEL, MLEV,LFINE '',3I5)')level,mlev,lfine
!
  23000 IF(.NOT. (level <= lfine))GO TO 23001
  mptr    = lstart(level)
  23002     IF(.NOT. (mptr /= 0))GO TO 23003
  CALL outmsh(mptr)
!c               if (dmpc) call dmpcnst(mptr)
  mptr = node(10, mptr)
  GO TO 23002
  23003     CONTINUE
  level = level + 1
  GO TO 23000
  23001 CONTINUE
!
  RETURN
END SUBROUTINE outtre
!
! --------------------------------------------------------------------
!

SUBROUTINE outmsh(mptr)
  INCLUDE 'nodal.inc'
!
!  this routine prints out the grid information
!
  IF( (node(4,mptr) < 1) ) RETURN
  time = rnode(20,mptr)/60.
  WRITE(6,1)mptr,node(4,mptr),time
  1     FORMAT(/,1X,' ******** GRID ',i3,'  LEVEL ',i3,                 &
               ' AT ',f10.3,' MINUTES ')
  WRITE(6,2) (rnode(i,mptr),i=1,8),rnode(29,mptr),rnode(30,mptr)
  2     FORMAT(/,10X,'corners        x              y   ',              &
            /,13X,'1',3X,2(2X,e13.6),/,13X,'2',3X,2(2X,e13.6),          &
            /,13X,'3',3X,2(2X,e13.6),/,13X,'4',3X,2(2X,e13.6),          &
            //,' z at bottom = ',e13.6,' z at top = ',e13.6,/)
  WRITE(6,3) (rnode(i,mptr),i=9,11),rnode(31,mptr)
  3     FORMAT(3X,'  dx =  ',e12.5,'   dy =  ',e12.5,                   &
               /,3X,'  dt =  ',e12.5,'   dz =  ',e12.5)
  WRITE(6,4)rnode(21,mptr),rnode(22,mptr),                              &
            (rnode(i,mptr),i=12,19),                                    &
            (rnode(i,mptr),i=23,28)
  4     FORMAT(3X,'  cos = ',e12.5,'   sin = ',e12.5,                   &
               //,'  misc grid description information  ',              &
               /,4(/,2X,4(1X,e12.5)) )
  WRITE(6,5) node(5,mptr),node(6,mptr),node(14,mptr),                   &
             node(7,mptr),node(8,mptr),node(16,mptr)
  5     FORMAT(/,'  NODAL INFORMATION ',/,                              &
                 1X,'   nx = ',i4,'  ny = ',i4,'  nz = ',i4,/,          &
                 1X,'   bigstep file = ',i3,' smlstp file = ',i3,       &
                    ' constants file = ',i3)
  WRITE(6,6) node(10,mptr),node(12,mptr),node(13,mptr)
  6     FORMAT(1X,'   next grid on level = ',i3,                        &
                  ' previous grid on level = ',i3,/,                    &
                1X,'   number of timesteps taken = ',i4)
  RETURN
END SUBROUTINE outmsh
!
! --------------------------------------------------------------------
!

SUBROUTINE outch( a,tmp,m,n,kk )
  DIMENSION a(m,n,kk),tmp(1)
!
!  plot out a few levels
!
  kmid = kk/2
  jmid = n/2
  imid = m/2
!  call parray(a(1,1,kmid),m,n)
  DO k=1,kk
    DO j=1,n
      tmp(j+(k-1)*n) = a(imid,j,k)
    END DO
  END DO
!  call parray(tmp,n,kk)
  DO k=1,kk
    DO i=1,m
      tmp(i+(k-1)*m) = a(i,jmid,k)
    END DO
  END DO
!  call parray(tmp,m,kk)
  RETURN
END SUBROUTINE outch
!
! --------------------------------------------------------------------
!

SUBROUTINE dmpdsc
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
!
!  print the generic grid description information if wanted
!
  IF( .NOT. verbose2 ) RETURN
!
  WRITE(6,'(''  GENERIC GRID INFORMATION '')')
  WRITE(6,1) nsint,nsreal,nx1d,ny1d,nz1d,                               &
             nxy2d,nxz2d,nyz2d,nxyz3d
!
  IF (nxy2d > 0) THEN
    DO i=1,nxy2d
      WRITE(6,2) i
      WRITE(6,3) stgxy(1,i),stgxy(2,i)
      WRITE(6,4) idmxy(1,i),idmxy(2,i)
      WRITE(6,5) ipkxy(i)
      WRITE(6,6) iupxy(1,i),iupxy(2,i)
      WRITE(6,7) ibdxy(1,i),ibdxy(2,i),ibdxy(3,i)
    END DO
  END IF
!
  IF (nxz2d > 0) THEN
    DO i=1,nxz2d
      WRITE(6,12) i
      WRITE(6,3) stgxz(1,i),stgxz(2,i)
      WRITE(6,4) idmxz(1,i),idmxz(2,i)
      WRITE(6,5) ipkxz(i)
      WRITE(6,6) iupxz(1,i),iupxz(2,i)
      WRITE(6,7) ibdxz(1,i),ibdxz(2,i)
    END DO
  END IF
!
  IF (nyz2d > 0) THEN
    DO i=1,nyz2d
      WRITE(6,22) i
      WRITE(6,3) stgyz(1,i),stgyz(2,i)
      WRITE(6,4) idmyz(1,i),idmyz(2,i)
      WRITE(6,5) ipkyz(i)
      WRITE(6,6) iupyz(1,i),iupyz(2,i)
      WRITE(6,7) ibdyz(1,i),ibdyz(2,i)
    END DO
  END IF
!
  IF (nxyz3d > 0) THEN
    DO i=1,nxyz3d
      WRITE(6,32) i
      WRITE(6,33) stgxyz(1,i),stgxyz(2,i),stgxyz(3,i)
      WRITE(6,34) idmxyz(1,i),idmxyz(2,i),idmxyz(3,i)
      WRITE(6,5)  ipkxyz(i)
      WRITE(6,6)  iupxyz(1,i),iupxyz(2,i)
      WRITE(6,7)  ibdxyz(1,i),ibdxyz(2,i),ibdxyz(3,i)
    END DO
  END IF
!
  1   FORMAT( /,2X,' nsint  = ',i5,                                     &
              /,2X,' nsreal = ',i5,                                     &
              /,2X,' nx1d   = ',i5,                                     &
              /,2X,' ny1d   = ',i5,                                     &
              /,2X,' nz1d   = ',i5,                                     &
              /,2X,' nxy2d  = ',i5,                                     &
              /,2X,' nxz2d  = ',i5,                                     &
              /,2X,' nyz2d  = ',i5,                                     &
              /,2X,' nxyz3d = ',i5  )
!
  2   FORMAT( /,'  data for 2-D xy variable ',i4 )
!
  3   FORMAT( '  x staggering ',f6.3,'  y staggering ',f6.3 )
  4   FORMAT( '  x max dimension ' ,i4,'  y max dimension ',i4 )
  5   FORMAT( '  pack variable - ',i2)
  6   FORMAT( '  update - ',i2,'  other vector ',i2)
  7   FORMAT( '  boundi - ',i2,'  other vector ',i2,                    &
              '  t-dt field - ',i2)
!
  12   FORMAT( /,'  data for 2-D xz variable ',i4 )
  22   FORMAT( /,'  data for 2-D yz variable ',i4 )
  32   FORMAT( /,'  data for 3-D xyz variable ',i4 )
  33   FORMAT( '  x staggering ',f6.3,'  y staggering ',f6.3,           &
               '  z staggering ',f6.3                         )
  34   FORMAT( '  x max dimension ' ,i4,'  y max dimension ',i4,        &
               '  z max dimension ',i4 )
  RETURN
END SUBROUTINE dmpdsc
!

SUBROUTINE settmp
  INCLUDE 'agrialloc.inc'
  INCLUDE 'manage.inc'
  INCLUDE 'agricst.inc'
!
  ntemp = lfree(2,1)
  IF( verbose4 ) THEN
    WRITE(6,'(''  IN SETTMP '')')
    WRITE(6,'(''  NTEMP = '',I7)') ntemp
  END IF
  RETURN
END SUBROUTINE settmp
!
! --------------------------------------------------------------------
!

SUBROUTINE dmpstr
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
!
!  output storage locations for the grids if necessary
!
  IF( .NOT. verbose4 ) RETURN
!
  level = 1
  10    IF(level > lfine) GO TO 999
  mptr = lstart(level)
  20    IF(mptr == 0) GO TO 100
!
  WRITE(6,'(/,''  STORAGE LOCATIONS FOR GRID '',I4)') mptr
  WRITE(6,1001) ipint(mptr),ipreal(mptr),ips1d(mptr),ipx(mptr),         &
                ipy(mptr),ipz(mptr)
!
  IF(nxy2d > 0) THEN
    WRITE(6,'('' STORAGE LOCATIONS  FOR XY VARIABLES '')')
    DO i=1,nxy2d
      WRITE(6,1002) i,ipxy(i,mptr)
    END DO
  END IF
  IF(nxz2d > 0) THEN
    WRITE(6,'('' STORAGE LOCATIONS  FOR XZ VARIABLES '')')
    DO i=1,nxz2d
      WRITE(6,1002) i,ipxz(i,mptr)
    END DO
  END IF
  IF(nyz2d > 0) THEN
    WRITE(6,'('' STORAGE LOCATIONS  FOR YZ VARIABLES '')')
    DO i=1,nyz2d
      WRITE(6,1002) i,ipyz(i,mptr)
    END DO
  END IF
  IF(nxyz3d > 0) THEN
    WRITE(6,'('' STORAGE LOCATIONS  FOR XYZ VARIABLES '')')
    DO i=1,nxyz3d
      WRITE(6,1002) i,ipxyz(i,mptr)
    END DO
  END IF
  1001  FORMAT(2X,5(2X,i9))
  1002  FORMAT(2X,i5,2X,i9)
!
  mptr = node(10,mptr)
  GO TO 20
  100   level = level+1
  GO TO 10
  999   RETURN
END SUBROUTINE dmpstr
!

SUBROUTINE setrot(mlevel)
  INCLUDE 'nodal.inc'
!
!
!  SETROT:  SET THE ROTATION FIELDS IN RNODE FOR ALL GRIDS AT
!        THE SAME LEVEL AS MPTR.  THESE FIELDS ARE
!        COS AND SIN OF THE ANGLE OF ROTATION OF THE GRID.
!        (THESE ARE SET ELSEWHERE).
!        ALSO SET SLOPE/INTERCEPT FORM FOR ALL 4 LINES OF THE RECT.
!       SET THE MAX/MIN PROJECTION ON EACH EIGENVECTOR FIELD.
!
  mptr = mlevel
  23000 IF(.NOT. (mptr /= 0))GO TO 23001
  dist = SQRT( (rnode(7,mptr)-rnode(1,mptr))**2                         &
              +(rnode(8,mptr)-rnode(2,mptr))**2 )
  rnode( 21,mptr)=(rnode(7,mptr)-rnode(1,mptr))/dist
  rnode( 22,mptr)=(rnode(8,mptr)-rnode(2,mptr))/dist
!
  rnode( 12,mptr) =  rnode(21,mptr)
  rnode( 13,mptr) =  rnode(22,mptr)
  rnode( 14,mptr) = -rnode(22,mptr)
  rnode( 15,mptr) =  rnode(21,mptr)
!
  delx = rnode(3,mptr) - rnode(1,mptr)
  dely = rnode(4,mptr) - rnode(2,mptr)
  IF(.NOT. (ABS(delx) > .0001))GO TO 23002
  rnode(23,mptr) = dely / delx
  GO TO 23003
  23002     CONTINUE
  rnode(23,mptr) = 1000000000
  23003     CONTINUE
  IF(.NOT. (ABS(rnode(23,mptr)) > .0001))GO TO 23004
  rnode(24,mptr) = -1.0 / rnode(23,mptr)
  GO TO 23005
  23004     CONTINUE
  rnode(24,mptr) = - 1000000000
  23005     CONTINUE
  rnode(25,mptr) = rnode(2,mptr)-rnode(23,mptr)*rnode(1,mptr)
  rnode(26,mptr) = rnode(2,mptr)-rnode(24,mptr)*rnode(1,mptr)
  rnode(27,mptr) = rnode(6,mptr) - rnode(23,mptr)*rnode(5,mptr)
  rnode(28,mptr) = rnode(6,mptr) - rnode(24,mptr)*rnode(5,mptr)
!
  rnode( 16,mptr) = rnode(5,mptr)*rnode(12,mptr) +rnode(6,mptr)         &
                       *rnode(13,mptr)+.0001
  rnode( 17,mptr) = rnode(1,mptr)*rnode(12,mptr) +rnode(2,mptr)         &
                       *rnode(13,mptr)-.0001
  rnode( 18,mptr) = rnode(5,mptr)*rnode(14,mptr) +rnode(6,mptr)         &
                       *rnode(15,mptr)+.0001
  rnode( 19,mptr) = rnode(1,mptr)*rnode(14,mptr) +rnode(2,mptr)         &
                       *rnode(15,mptr)-.0001
  mptr = node(10, mptr)
  GO TO 23000
  23001 CONTINUE
!
  RETURN
END SUBROUTINE setrot
!

SUBROUTINE cpyfld( mptr,ifldin,ifldout,idim )
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
!
! this routine copies the data in xy or xyz field
! ifldin onto ifldout and replaces ifldout.
!
  nx = node(5,mptr)
  ny = node(6,mptr)
  nz = node(14,mptr)
!
  IF (idim == 2) THEN
    inptrs = igtnxy( mptr,  ifldin, 1 )
    inptrr = igtnxy( mptr, ifldout, 1 )
!
    DO ij=1,nx*ny
      a(inptrr+ij-1) = a(inptrs+ij-1)
    END DO
!
    CALL retnxy( mptr, ifldout, 1, inptrr, .true. )
!
  ELSE IF (idim == 3) THEN
!
    inptrs = igtxyz( mptr,  ifldin, 1 )
    inptrr = igtxyz( mptr, ifldout, 1 )
    DO ijk=1,nx*ny*nz
      a(inptrr+ijk-1) = a(inptrs+ijk-1)
    END DO
    CALL retxyz( mptr, ifldout, 1, inptrr, .true. )
!
  ELSE
    WRITE(6,'(''  ERROR, NOT 2-D OR 3-D IN CPYFLD '',                   &
    &             4I9)') mptr,ifldin,ifldout,idim
  END IF
  RETURN
END SUBROUTINE cpyfld

SUBROUTINE addgrd( mptr )

  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
  CHARACTER (LEN=256) :: filenm
!
!  this routine takes a grid stored at unit iunit
!  and adds it to the end of the grid storage
!
  CALL resett
  CALL setstr( mptr )
  CALL settmp
  in = ipint( mptr )
!
!  get the proper file
!
  CALL mkflnm( 'tmpdmp',' ','.z',0.,mptr,filenm,length )
  CALL fexist( filenm(1:length),length,.true.,.false. )
  iunit = igtunit( dum )
  OPEN( UNIT=iunit,FILE=filenm(1:length),FORM='unformatted',            &
        STATUS='old'                                        )
  IF( verbose3 ) WRITE(6,*) ' read grid ',mptr,' from tmp file ',       &
                filenm(1:length),', io through unit ',iunit
  REWIND( iunit )
!
! Yuhe: when EXBC turned on, we must add the EXBC arrays for grid 1.
!
  IF ( lexbc == 1 .AND. mptr == 1 ) THEN
    nwords = ipexbc(nexbc3d+1,mptr) - ipint(mptr)
  ELSE
    nwords = ipxyz(nxyz3d+1,mptr) - ipint(mptr)
  END IF

  IF( verbose3 ) WRITE(6,*) ' calling rdngrd ',mptr,iunit,nwords
  CALL rdngrd( iunit,a(in),nwords )
!
!  close and delete file, return unit number
!
  CLOSE( UNIT=iunit,STATUS='keep' )
  CALL retnunit( iunit )
!
! we're finished here
!
  RETURN
END SUBROUTINE addgrd
!

SUBROUTINE dmpgrd( mptr,saved )

  LOGICAL :: saved
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
  CHARACTER (LEN=256) :: filenm

!
!  this routine removes grid storage from the incore
!  array and, if desired, writes it out to file iunit
!  which, at present, will be the grid number plus 39
!
  IF ( saved ) THEN
!
! dump data to disk file
!
    CALL mkflnm( 'tmpdmp',' ','.z',0.,mptr,filenm,length )
    iunit = igtunit( dum )
    CALL fexist( filenm(1:length),length,.false.,.false. )
    OPEN( UNIT=iunit,FILE=filenm(1:length),FORM='unformatted',          &
          STATUS='new'                                        )
    IF( verbose3 ) WRITE(6,*) ' dump grid ',mptr,' to tmp file ',       &
                 filenm(1:length),', io through unit ',iunit
    REWIND( iunit )
!
! Yuhe: when EXBC turned on, we must add the EXBC arrays for grid 1.
!
    IF ( lexbc == 1 .AND. mptr == 1 ) THEN
      nwords = ipexbc(nexbc3d+1,mptr) - ipint(mptr)
    ELSE
      nwords = ipxyz(nxyz3d+1,mptr) - ipint(mptr)
    END IF

    IF( verbose3 ) WRITE(6,*) ' calling wrtgrd ',mptr,iunit,nwords
    CALL wrtgrd( iunit, a(ipint(mptr)), nwords )
!
!  close and delete file, return unit number
!
    CLOSE( UNIT=iunit,STATUS='keep' )
    CALL retnunit( iunit )
!
  END IF
!
! free up this space
!
  IF( verbose4 ) WRITE(6,*) ' reclaiming space for grid ',mptr
  CALL resett
  CALL reclam( ipint(mptr),nwords )
  CALL settmp
!
! we're finished here
!
  RETURN
END SUBROUTINE dmpgrd
!

SUBROUTINE rdngrd( iunit,a,nwords )
  DIMENSION a(nwords)
  READ(iunit) a
  RETURN
END SUBROUTINE rdngrd
!

SUBROUTINE wrtgrd( iunit,a,nwords )
  DIMENSION a(nwords)
  WRITE(iunit) a
  RETURN
END SUBROUTINE wrtgrd
!

SUBROUTINE dmplvl( level )
  INCLUDE 'nodal.inc'
  INCLUDE 'agricst.inc'
!
!  this dumps out all the grids at level 'level'
!  into individual files using the dmpgrd subroutine
!
!  grids in lstart are dumped, starting with the last one.
!  it is assumed that the present level is the highest
!  existing in memory
!
  mptr = lback(level)
  IF(mptr == 0) THEN
    WRITE(6,*) ' error in dmplvl, no grids on level ',level
    STOP
  END IF
  5     IF( mptr == 0) GO TO 99
  IF( verbose3 ) WRITE(6,*) ' dumping grid - level ',mptr,level
  CALL dmpgrd( mptr,.true. )
  mptr = node(12,mptr)
  GO TO 5
!
  99    RETURN
END SUBROUTINE dmplvl
!

SUBROUTINE exchng( level )
  INCLUDE 'nodal.inc'
  INCLUDE 'agricst.inc'
!
!  this subroutine calls the boundary condition exchange
!  for overlapping fine grids
!
!  this version doesn't do an overlap check first.
!  later version should, we just need to put in overlap
!  test to do so
!
  cpu0 = f_cputime()
  mptr = lstart( level )
  5     IF(mptr == 0) GO TO 999
!
  mptre = node(10,mptr)
  10    IF(mptre == 0) GO TO 20
!
!  do boundary condition exchange between mptr and mptre
!
!  here's where the overlap test will be
!
  WRITE(6,*) ' b.c. exchange between ',mptr,mptre
  CALL exchbc( mptr,mptre )
  mptre = node(10,mptre)
  GO TO 10
!
  20    CONTINUE
!
!  we've gone through the list for grid mptr, now go
!  to the next grid and go through the rest of the list
!
  mptr = node(10,mptr)
  GO TO 5
!
  999   CONTINUE
  cpu_bndexch = cpu_bndexch + f_cputime() - cpu0
!
!  we're finished here
!
  RETURN
END SUBROUTINE exchng
!

SUBROUTINE dmpall( time )
  INCLUDE 'nodal.inc'
  INCLUDE 'agricst.inc'
!
!  this dumps out all the grids for a restart
!  into individual files using the dumpgrid subroutine.
!  this routine does not clear internal storage, thus integrations
!  can be continued after this dump, as opposed to the
!  dumps performed by dmplvl, which does clear the internal storage
!
!  grids in lback are dumped, starting with the last one.
!
  level = 1
  mptr = lback(level)
!
  IF(mptr == 0) THEN
    WRITE(6,*) ' error in dmpall, no grids on level 1 '
    STOP
  END IF
!
  5    IF( mptr == 0) GO TO 10
!
  IF( verbose3 ) WRITE(6,*) ' dumping grid - level ',mptr,level
  CALL dumpgrid( mptr,time )
  mptr = node(12,mptr)
  GO TO 5
  10    level = level + 1
  mptr = lback(level)
  IF(level <= lfine) GO TO 5
!
  99    RETURN
END SUBROUTINE dmpall
!

SUBROUTINE readall( time )
  INCLUDE 'nodal.inc'
  INCLUDE 'agricst.inc'
!
!  this reads in all the grids for a restart
!  from individual files using the readgrid subroutine.
!  this routine assumes that the internal storage is
!  already configured, thus no storage management needs to be
!  done (as opposed to reading in a new grid using addgrd,
!  which does allocate new mwmory for the grid)
!
  level = 1
  mptr = lback(level)
!
  IF(mptr == 0) THEN
    WRITE(6,*) ' error in readall, no grids on level 1 '
    STOP
  END IF
!
  5   IF( mptr == 0) GO TO 10
!
  IF( verbose3 ) WRITE(6,*) ' reading grid - level ',mptr,level
  CALL readgrid( mptr,time )
  mptr = node(12,mptr)
  GO TO 5
  10    level = level + 1
  mptr = lback(level)
  IF(level <= lfine) GO TO 5
!
  99    RETURN
END SUBROUTINE readall
!
!

SUBROUTINE dumpgrid( mptr,time )

  IMPLICIT NONE

  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
  INCLUDE 'globcst.inc'

  CHARACTER (LEN=256) :: filenm
  INTEGER :: length

  INTEGER :: mptr
  INTEGER :: nwords
  INTEGER :: iunit, igtunit
  INTEGER :: machsti

  REAL :: dum
  REAL :: time
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  this routine dumps out grid data for an individual
!  grid for later restart.
!
  CALL mkflnm( runname(1:lfnkey),' ','.t',time,mptr,filenm,length )

  iunit = igtunit( dum )
  WRITE(6,*) ' restart dump for run ',runname(1:lfnkey),                &
             ', grid ',mptr,' at time ',time

  IF( verbose3 ) WRITE(6,*) ' dumping to file ',filenm(1:length),       &
                ' through unit ',iunit
  CALL fexist( filenm(1:length),length,.false.,.true. )

  OPEN( UNIT=iunit,FILE=filenm(1:length),FORM='unformatted',            &
        STATUS='new' )
  REWIND( iunit )
  WRITE(iunit) runname,filenm,machst,lfnkey,mptr,time

!
! Yuhe: when EXBC turned on, we must add the EXBC arrays for grid 1.
!
  IF ( lexbc == 1 .AND. mptr == 1 ) THEN
    nwords = ipexbc(nexbc3d+1,mptr) - ipint(mptr)
  ELSE
    nwords = ipxyz(nxyz3d+1,mptr) - ipint(mptr)
  END IF

  IF( verbose3 ) WRITE(6,*) ' calling wrtgrd ',mptr,iunit,nwords
  CALL wrtgrd( iunit, a(ipint(mptr)), nwords )
  CLOSE( UNIT=iunit,STATUS='keep' )

  CALL retnunit( iunit )
  IF( verbose3 ) WRITE(6,*) ' closed ',filenm(1:length),                &
                ' on unit ',iunit
!
! we're finished dumping this grid data
!
  RETURN
END SUBROUTINE dumpgrid
!

SUBROUTINE readgrid( mptr,time )

!  implicit none

  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
!  include 'globcst.inc'

  CHARACTER (LEN=256) :: filenm

!
!  this routine reads in grid data for an individual
!  grid for restart.
!
  CALL mkflnm( runold(1:nmlntho),' ','.t',time,mptr,filenm,length )
  iunit = igtunit( dum )
  WRITE(6,*) ' restart read from run ',runold(1:nmlntho),               &
             ', grid ',mptr,' at time ',time
  IF( verbose3 ) WRITE(6,*) ' reading from file ',filenm(1:length),     &
                ' through unit ',iunit
  CALL fexist( filenm(1:length),length,.true.,.false. )

  OPEN( UNIT=iunit,FILE=filenm(1:length),FORM='unformatted',            &
        STATUS='old' )
  REWIND( iunit )
  READ( iunit ) runnami,filenmi,machsti,nmlnthi,mptri,timei
!
!  here we assume everything is ok, we could check the
!  runname and filenames etc... that come from the file
!  if that is deemed necessary
!
!
! Yuhe: when EXBC turned on, we must add the EXBC arrays for grid 1.
!
  IF ( lexbc == 1 .AND. mptr == 1 ) THEN
    nwords = ipexbc(nexbc3d+1,mptr) - ipint(mptr)
  ELSE
    nwords = ipxyz(nxyz3d+1,mptr) - ipint(mptr)
  END IF

  IF( verbose3 ) WRITE(6,*) ' calling rdngrd ',mptr,iunit,nwords
  CALL rdngrd( iunit, a(ipint(mptr)), nwords )

  CLOSE( UNIT=iunit,STATUS='keep' )
  CALL retnunit( iunit )
  IF( verbose3 ) WRITE(6,*) ' closed ',filenm(1:length),                &
                ' on unit ',iunit
!
! we're finished reading this grid data
!
  RETURN
END SUBROUTINE readgrid
!
!

SUBROUTINE iosetup
  COMMON / unitnmb / iunitn(10),nleft,ntot
  DO i=1,10
    iunitn(i) = 40+i
  END DO
  nleft = 10
  ntot = 10

  WRITE(6,'(/2(/1x,a)/)')                                               &
      'Fortran I/O units 40-50 reserved for AGR interface.',            &
      'Do not use them for other purposes.'

  RETURN
END SUBROUTINE iosetup
!
!

  INTEGER FUNCTION igtunit( dum )
  COMMON / unitnmb / iunitn(10),nleft,ntot
  INCLUDE 'agricst.inc'
!
  IF(nleft == 0) THEN
    WRITE(6,*) ' no more units left for io in iunitn list '
    WRITE(6,*) ' error stop '
  END IF
!
  igtunit = iunitn( nleft )
  IF( verbose3 ) WRITE(6,*) ' igtunit getting unit= ',igtunit
  nleft = nleft-1
  RETURN
  END FUNCTION igtunit
!
!

SUBROUTINE retnunit( iunit )
  COMMON / unitnmb / iunitn(10),nleft,ntot
!
  nleft = nleft+1
  IF(nleft > ntot) THEN
    WRITE(6,*) ' something is screwed, there are '
    WRITE(6,*) ' more units in iunitn than when we started!! '
    WRITE(6,*) ' returning unit ',iunit
    WRITE(6,*) ' ntot, nleft = ',ntot,nleft
    DO i=1,ntot
      WRITE(6,2) i,iunitn(i)
    END DO
  END IF
  2    FORMAT(2X,i5,1X,i5)
!
  iunitn(nleft) = iunit
  RETURN
END SUBROUTINE retnunit
!
!-----------------------------------------------------------------------
!

SUBROUTINE mkflnm( instr1,instr2,instr3,                                &
           time,grid,outstr,length )

  IMPLICIT NONE
  REAL :: time
  INTEGER :: i, j, k, length, grid, tlength,glength
  CHARACTER (LEN=*) :: instr1, instr2, instr3, outstr
  CHARACTER (LEN=80) :: tmp, chmake
  EXTERNAL chmake
  DATA tlength,glength/ 5, 2 /

!-----------------------------------------------------------------------
! Find lengths of the strings
!-----------------------------------------------------------------------
  i = INDEX( instr1,' ' ) - 1
  IF(i == -1) i = LEN( instr1 )
  j = INDEX( instr2,' ' ) - 1
  IF(j == -1) j = LEN( instr2 )
  k = INDEX( instr3,' ' ) - 1
  IF(k == -1) k = LEN( instr3 )
  length = 0
!
  IF( i > 0 ) THEN
    outstr(1:i) = instr1(1:i)
    length = i
  END IF
!
  IF( j > 0 ) THEN
    outstr(length+1:length+j) = instr2(1:j)
    length = length+j
  END IF
!
!-----------------------------------------------------------------------
! if time is < 0, then all we want is this string
!-----------------------------------------------------------------------
!
  IF( time < 0 ) RETURN
!
!-----------------------------------------------------------------------
! If TIME > 0 add time to the string
!-----------------------------------------------------------------------
!
  tmp = chmake( IFIX(time),tlength )
!
  IF(k > 0) THEN
    outstr(i+j+1:) = instr3(1:k)//tmp(1:tlength)
  ELSE
    outstr(i+j+1:) = tmp(1:tlength)
  END IF
!
  length = length + tlength + k
!
!-----------------------------------------------------------------------
! Create the grid no string if grid number greater than zero
!-----------------------------------------------------------------------
!
  IF( grid > 0 ) THEN
    tmp = chmake( grid,glength )
!      outstr(length+1:) = '.g'//tmp(1:glength)
!      length = length + glength + 2
    outstr(length+1:) = '.'//tmp(1:glength)
    length = length + glength + 1
  END IF
!
!  we're finished
!
  WRITE(6,'(a,a)') 'The file name created is ',outstr(1:length)

  999    RETURN
END SUBROUTINE mkflnm
!
!

  FUNCTION chmake( numb,length )

  IMPLICIT NONE
  INTEGER :: numb,length,numbc,i
  CHARACTER (LEN=*) :: chmake
  CHARACTER (LEN=16) :: tmp

!
!  check for degenerate cases
!
  IF(length <= 0)  THEN
    WRITE(6,*) ' length is <= zero in chmake ',length
    WRITE(6,*) ' error exit '
    STOP
  END IF

!
  IF(numb == 0) THEN
!
!  make the string
!
    DO i=1,length
      chmake(i:i) = '0'
    END DO
    RETURN
  END IF
!
  numbc = INT(ALOG10(ABS(FLOAT(numb)))) + 1
  IF( numb < 0 ) numbc = numbc+1
!
  IF( numbc > length ) THEN
    WRITE(6,*) ' number too large in chmake ',                          &
               ' for given length ',length,numb,numbc
    WRITE(6,*) ' error exit '
    STOP
  END IF
!
  IF( numbc > 16 ) THEN
    WRITE(6,*) ' number has more than 16 chars ',                       &
               ' in chmake, you must be kidding!!! '
    WRITE(6,*) ' error exit '
    STOP
  END IF
!
!  make the string
!
  DO i=1,length
    chmake(i:i) = '0'
  END DO
!
  WRITE(tmp,2) numb
  1      FORMAT(i1)
  2      FORMAT(i16)
  chmake(length-numbc+1:length) = tmp(16-numbc+1:16)
!
!  we're finished here
!
  RETURN
  END FUNCTION chmake
!
!

SUBROUTINE cleanf( mptrs )
  INCLUDE 'nodal.inc'
  INCLUDE 'agricst.inc'
  CHARACTER (LEN=256) :: filenm
  LOGICAL :: existf
!
!  this routine closes the files for the temp grid
!  data dumps used in regrid, mptr is the first grid
!  on  a level to be dumped
!
  iunit = igtunit( dum )
  mptr = mptrs
  10   IF( mptr == 0) GO TO 99
  CALL mkflnm( 'tmpdmp',' ','.z',0.,mptr,filenm,length )
  INQUIRE(FILE=filenm(1:length),EXIST=existf)
!
  IF(existf) THEN
    OPEN( UNIT=iunit,FILE=filenm(1:length),                             &
          FORM='unformatted',STATUS='old'   )
    CLOSE(UNIT=iunit,STATUS='delete')
  ELSE
    WRITE(6,*) ' in cleanf, expected to find file ',                    &
                filenm(1:length),', yet does not exist '
  END IF
  mptr = node(10,mptr)
  GO TO 10
  99   CONTINUE
  CALL retnunit( iunit )
  RETURN
END SUBROUTINE cleanf
!

SUBROUTINE fexist( filenm,length,exst,incr )
  CHARACTER (LEN=*) :: filenm
  INTEGER :: length
  LOGICAL :: exst,incr,existf
!
!  this subroutine checks on the existence of the file
!  filenm(1:length).
!
!  exst: true if the file should exist, false if it shouldn't
!     exist
!  incr: if true and exst is false but the file exists, then
!     return a new name for a file that doesn't exist.
!     here we just patch on a version number to the file
!
!  if the file should exist but doesn't, error exit
!   if the file shouldn't exist but does and no increment
!   is wanted, then we error exit
!   later we'll build in an erase option taht will allow
!   removal of existing files that we do not want.
!
!  first we check for file that should exist
!
  IF( exst ) THEN
    INQUIRE( FILE=filenm(1:length),EXIST=existf )
    IF(.NOT.existf ) THEN
      WRITE(6,*) ' exected to find file ',filenm(1:length)
      WRITE(6,*) ' file does not exist, error exit '
      STOP
    ELSE
      RETURN
    END IF
  END IF
!
!  now for file that shouldn't exist
!
  inc = 0
  10   INQUIRE( FILE=filenm(1:length),EXIST=existf )
!
  IF (.NOT. existf ) RETURN
!
  IF (.NOT. incr ) THEN
    WRITE(6,*) ' cannot overwrite pre-existing file ',                  &
                filenm(1:length),', error exit '
    STOP
  END IF
!
!  here we put a number onto the end of the file
!
  inc = inc+1
  IF(inc == 1) THEN
    WRITE(filenm(length+1:length+2),15) inc
    length = length+2
  ELSE IF(inc < 10) THEN
    WRITE(filenm(length-1:length  ),15) inc
  ELSE IF(inc > 10) THEN
    WRITE(6,*) '  to many version numbers of file ',                    &
               filenm(1:length-2),', error exit '
    STOP
  END IF
  GO TO 10
!
  999  RETURN
  15   FORMAT('.',i1)
END SUBROUTINE fexist

SUBROUTINE filzero(a,n)
  REAL :: a(n)

  DO i=1,n
    a(i) = 0.0
  END DO

  RETURN
END SUBROUTINE filzero
!
!-----------------------------------------------------------------------
!
! Created by Louis Wicker, November 1992
! Writes out the cpu usage for the model
!
! Modified for ARPS by Ming Xue, 12/10/1992.
!
!-----------------------------------------------------------------------
!

SUBROUTINE prtcpu(flag,UNIT)

  INCLUDE 'agricpu.inc'
  INTEGER :: flag, UNIT

!-----------------------------------------------------------------------
!  Flag = 0, initialize values to zero
!-----------------------------------------------------------------------

  IF( flag == 0 ) THEN

    cpu_main       = 0.0
    cpu_advect     = 0.0
    cpu_smlstep    = 0.0
    cpu_mixing     = 0.0
    cpu_filter     = 0.0
    cpu_cloud      = 0.0
    cpu_boundary   = 0.0
    cpu_copy       = 0.0
    cpu_init0      = 0.0
    cpu_usrout     = 0.0
    cpu_regrid     = 0.0
    cpu_bndexch    = 0.0
    cpu_bndcint    = 0.0
    cpu_update     = 0.0
    cpu_pack       = 0.0
    cpu_unpack     = 0.0
    cpu_simulation = 0.0

!-----------------------------------------------------------------------
!  Flag=1, write out cpu information
!-----------------------------------------------------------------------

  ELSE

    cpu_solver = cpu_advect+cpu_smlstep+cpu_mixing+cpu_cloud            &
               + cpu_boundary+cpu_copy+cpu_main+cpu_filter

    cpu_interface = cpu_simulation - cpu_solver - cpu_usrout            &
                  - cpu_init0

    WRITE(UNIT,999)
    WRITE(UNIT,101) 'ADAPTIVE MESH CLOUD MODEL CPU STATISTICS'
    WRITE(UNIT,999)

    WRITE(UNIT,101) 'Total CPU for SIMULATION  is ',cpu_simulation
    WRITE(UNIT,999)

    WRITE(UNIT,101) 'Total CPU for INTERFACE   is ',cpu_interface
    WRITE(UNIT,101) 'Total CPU for REGRID      is ',cpu_regrid
    WRITE(UNIT,101) 'Total CPU for EXCHBC      is ',cpu_bndexch
    WRITE(UNIT,101) 'Total CPU for UPDBC       is ',cpu_bndcint
    WRITE(UNIT,101) 'Total CPU for UPDATE      is ',cpu_update
    WRITE(UNIT,101) 'Total CPU for PACKING     is ',cpu_pack
    WRITE(UNIT,101) 'Total CPU for UNPACKING   is ',cpu_unpack
    WRITE(UNIT,101) 'Total CPU for INIT        is ',cpu_init0
    WRITE(UNIT,101) 'Total CPU for USROUT      is ',cpu_usrout
    WRITE(UNIT,999)

    WRITE(UNIT,101) 'Total CPU for SOLVER      is ',cpu_solver
    WRITE(UNIT,101) 'Total CPU for CLOUD3D     is ',cpu_main
    WRITE(UNIT,101) 'Total CPU for SUB SMLSTEP is ',cpu_smlstep
    WRITE(UNIT,101) 'Total CPU for SUB ADVECT  is ',cpu_advect
    WRITE(UNIT,101) 'Total CPU for SUB CLOUD   is ',cpu_cloud
    WRITE(UNIT,101) 'Total CPU for SUB MIXING  is ',cpu_mixing
    WRITE(UNIT,101) 'Total CPU for SUB FILTER  is ',cpu_filter
    WRITE(UNIT,101) 'Total CPU for SUB BOUNDS  is ',cpu_boundary
    WRITE(UNIT,101) 'Total CPU for SUB COPY    is ',cpu_copy

    WRITE(UNIT,999)
    101     FORMAT(1X,a,f15.5)
    999     FORMAT(/79('-')/)

  END IF

  RETURN
END SUBROUTINE prtcpu


SUBROUTINE writint( icnst,ncnst , UNIT)

  INTEGER :: ncnst , UNIT
  INTEGER :: icnst(ncnst)

  DO i=1,ncnst
    WRITE(UNIT, '(1x,2i10)') i, icnst(i)
  END DO

  RETURN
END SUBROUTINE writint

SUBROUTINE readint( icnst,ncnst , UNIT)

  INTEGER :: ncnst , UNIT
  INTEGER :: icnst(ncnst)

  DO i=1,ncnst
    READ(UNIT, '(1x,2i10)') ii, icnst(i)
    WRITE(6, '(1x,2(a,i10))') 'Integer No. ',i,'=',icnst(i)
  END DO

  RETURN
END SUBROUTINE readint


