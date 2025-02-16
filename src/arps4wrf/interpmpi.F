!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE FIND_ARPSBDRWIDTH           ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE find_arpsbdrwidth(mm1s,mm1e,mm2s,mm2e,mu1s,mu1e,mu2s,mu2e,   &
                        mv1s,mv1e,mv2s,mv2e,pm1s,pm1e,pm2s,pm2e,        &
                        pu1s,pu1e,pu2s,pu2e, pv1s,pv1e,pv2s,pv2e,       & 
                        nx,ny,dx, dy,                                   &
                        mx,my,ux,uy,vx,vy,                              &
                        x,y,xs,ys,                                      &
                        lvldbg,bdrwidth,istatus)
  IMPLICIT NONE

!#######################################################################
  ! to ensure that 16 pt interpolation can be performed.

  INTEGER, INTENT(IN) :: mm1s,mm1e,mm2s,mm2e, mu1s,mu1e,mu2s,mu2e, mv1s,mv1e,mv2s,mv2e
  INTEGER, INTENT(IN) :: pm1s,pm1e,pm2s,pm2e, pu1s,pu1e,pu2s,pu2e, pv1s,pv1e,pv2s,pv2e
  REAL,    INTENT(IN) :: mx(mm1s:mm1e,mm2s:mm2e), my(mm1s:mm1e,mm2s:mm2e)
  REAL,    INTENT(IN) :: ux(mu1s:mu1e,mu2s:mu2e), uy(mu1s:mu1e,mu2s:mu2e)
  REAL,    INTENT(IN) :: vx(mv1s:mv1e,mv2s:mv2e), vy(mv1s:mv1e,mv2s:mv2e)

  INTEGER, INTENT(IN) :: nx, ny
  REAL,    INTENT(IN) :: xs(nx), ys(ny), x(nx), y(ny)
  REAL,    INTENT(IN) :: dx, dy  
  INTEGER, INTENT(IN) :: lvldbg
  
  INTEGER, INTENT(OUT) :: bdrwidth
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

  REAL    :: xmin,ymin,xmax,ymax

  INTEGER :: bdrxms, bdrxme, bdryms, bdryme
  INTEGER :: bdrxus, bdrxue, bdryus, bdryue
  INTEGER :: bdrxvs, bdrxve, bdryvs, bdryve

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  bdrwidth = 0
#ifdef MPI_PARALLEL

!---------------------- M grid -----------------------------------------

  xmin = MINVAL( mx(pm1s:pm1e,pm2s:pm2e) )
  xmax = MAXVAL( mx(pm1s:pm1e,pm2s:pm2e) )

  ymin = MINVAL( my(pm1s:pm1e,pm2s:pm2e) )
  ymax = MAXVAL( my(pm1s:pm1e,pm2s:pm2e) )

  IF (lvldbg > 2) THEN
    WRITE(6,'(1x,a,4F15.2)') 'Mass Points ARPS range: ',xs(2),xs(nx-2),ys(2),ys(ny-2)
    WRITE(6,'(1x,a,4F15.2)') 'Mass Points WRF  range: ',xmin,xmax,ymin,ymax
  END IF
  
  bdrxms = 0
  bdrxme = 0
  bdryms = 0
  bdryme = 0

  bdrxus = 0
  bdrxue = 0
  bdryus = 0
  bdryue = 0

  bdrxvs = 0
  bdrxve = 0
  bdryvs = 0
  bdryve = 0

  IF( xmin < xs(2) ) THEN
    bdrxms = NINT( (xs(2)-xmin)/dx+0.5 )
  END IF
  IF( xmax > xs(nx-2) ) THEN
    bdrxme = NINT( (xmax-xs(nx-2))/dx+0.5 )
  END IF

  IF( ymin < ys(2) ) THEN
    bdryms = NINT( (ys(2)-ymin)/dy+0.5 )
  END IF

  IF( ymax > ys(ny-2) ) THEN
    bdryme = NINT( (ymax-ys(ny-2))/dy+0.5 )
  END IF

  bdrwidth = MAX( bdrxms, bdrxme, bdryms, bdryme ) 

  IF (lvldbg > 1) THEN
   WRITE(6,'(1x,a,I2)') 'Mass points message boundary = ', bdrwidth
  END IF

!---------------------- U grid -----------------------------------------

  xmin = MINVAL( ux(pu1s:pu1e,pu2s:pu2e) )
  xmax = MAXVAL( ux(pu1s:pu1e,pu2s:pu2e) )

  ymin = MINVAL( uy(pu1s:pu1e,pu2s:pu2e) )
  ymax = MAXVAL( uy(pu1s:pu1e,pu2s:pu2e) )
  
  IF (lvldbg > 2) THEN
    WRITE(6,'(1x,a,4F15.2)') 'U Points ARPS range: ',x(2),x(nx-1),ys(2),ys(ny-2)
    WRITE(6,'(1x,a,4F15.2)') 'U Points WRF  range: ',xmin,xmax,ymin,ymax
  END IF
  
  IF( xmin < x(2) ) THEN
    bdrxus = NINT( (x(2)-xmin)/dx+0.5 )
  END IF

  IF( xmax > x(nx-1) ) THEN
    bdrxue = NINT( (xmax-x(nx-1))/dx+0.5 )
  END IF

  IF( ymin < ys(2) ) THEN
    bdryus = NINT( (ys(2)-ymin)/dy+0.5 )
  END IF

  IF( ymax > ys(ny-2) ) THEN
    bdryue = NINT( (ymax-ys(ny-2))/dy+0.5 )
  END IF

  bdrwidth = MAX( bdrwidth, bdrxus, bdrxue, bdryus, bdryue ) 

  IF (lvldbg > 1) THEN
   WRITE(6,'(1x,a,I2)') 'U points message boundary = ', bdrwidth
  END IF

!---------------------- V grid -----------------------------------------

  xmin = MINVAL( vx(pv1s:pv1e,pv2s:pv2e) )
  xmax = MAXVAL( vx(pv1s:pv1e,pv2s:pv2e) )

  ymin = MINVAL( vy(pv1s:pv1e,pv2s:pv2e) )
  ymax = MAXVAL( vy(pv1s:pv1e,pv2s:pv2e) )
  
  IF (lvldbg > 2) THEN
    WRITE(6,'(1x,a,4F15.2)') 'V Points ARPS range: ',xs(2),xs(nx-2),y(2),y(ny-1)
    WRITE(6,'(1x,a,4F15.2)') 'V Points WRF  range: ',xmin,xmax,ymin,ymax
  END IF
  
  IF( xmin < xs(2) ) THEN
    bdrxvs = NINT( (xs(2)-xmin)/dx+0.5 )
  END IF

  IF( xmax > xs(nx-2) ) THEN
    bdrxve = NINT( (xmax-xs(nx-2))/dx+0.5 )
  END IF

  IF( ymin < y(2) ) THEN
    bdryvs = NINT( (y(2)-ymin)/dy+0.5 )
  END IF

  IF( ymax > y(ny-1) ) THEN
    bdryve = NINT( (ymax-y(ny-1))/dy+0.5 )
  END IF

  bdrwidth = MAX( bdrwidth, bdrxvs, bdrxve, bdryvs, bdryve ) 

  IF (lvldbg > 1) THEN
   WRITE(6,'(1x,a,I2)') 'V points message boundary = ', bdrwidth
  END IF

  CALL mpmaxi(bdrwidth)

!  bdrwidth = MAX(bdrwidth-1,0)   ! ARPS arrays contain 1 fake zone
!                                 ! on each side.
#endif
  RETURN
END SUBROUTINE find_arpsbdrwidth
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SET_COORD_HALO              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE set_coord_halo(mxs,mxe,mys,mye,nx,ny,                        &
                      x,y,xs,ys,dx, dy,                                 &
                      xext,yext,xsext,ysext,istatus)

!#######################################################################
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: mxs, mxe, mys, mye, nx, ny
  REAL,    INTENT(IN)  :: x(nx),y(ny),xs(nx),ys(ny)
  REAL,    INTENT(IN)  :: dx,dy
  REAL,    INTENT(OUT) :: xext(mxs:mxe), yext(mys:mye)
  REAL,    INTENT(OUT) :: xsext(mxs:mxe),ysext(mys:mye)
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

#ifdef MPI_PARALLEL
  INCLUDE 'mpif.h'
  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
#endif

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_w for west boundary
                             ! mptag + tag_e for east boundary

  INTEGER :: bdrwidth
  INTEGER :: i,j,n

  REAL, ALLOCATABLE :: tem1(:), tem2(:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  bdrwidth = 1-mxs

  DO i = 1,nx
    xext(i)  = x(i)
    xsext(i) = xs(i)
  END DO

  DO j = 1,ny
    yext(j)  = y(j)
    ysext(j) = ys(j)
  END DO

  IF (loc_x == 1) THEN
    DO i = mxs,0
      n = 1-i
      xext(i)  =  x(1)-n*dx
      xsext(i) = xs(1)-n*dx
    END DO
  END IF

  IF (loc_x == nproc_x) THEN
    DO i = nx+1,mxe
      n = i-nx
      xext(i)  =  x(nx)+n*dx
      xsext(i) = xs(nx)+n*dx
    END DO
  END IF

  IF (loc_y == 1) THEN
    DO j = mys,0
      n = 1-j
      yext(j)  =  y(1)-n*dy
      ysext(j) = ys(1)-n*dy
    END DO
  END IF

  IF (loc_y == nproc_y) THEN
    DO j = ny+1,mye
      n = j-ny
      yext(j)  =  y(ny)+n*dy
      ysext(j) = ys(ny)+n*dy
    END DO
  END IF

#ifdef MPI_PARALLEL

  IF (bdrwidth <= 0) RETURN

  ALLOCATE(tem1(bdrwidth+1), STAT = istatus)
  ALLOCATE(tem2(bdrwidth+1), STAT = istatus)

  CALL inctag
  mptag = gentag

!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == nproc_x) THEN       ! last processor in a row
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  ! 
  ! receive from
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

!---------------------------- U-GRID -----------------------------------

  !
  ! Pack send buffer, send east slice to update west boundary of 
  ! the east neighbor
  !
  DO i=1,bdrwidth
    tem1(i) = x(nx-3-bdrwidth+i)
  END DO

  CALL mpi_sendrecv(tem1,bdrwidth,MPI_REAL,dest,  mptag+tag_w,       &
                    tem2,bdrwidth,MPI_REAL,source,mptag+tag_w,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update West boundary data
  !
  IF ( loc_x /= 1 ) THEN
    DO i=1,bdrwidth
      xext(i-bdrwidth) = tem2(i)
    END DO
  END IF

!---------------------------- MASS -----------------------------------

  !
  ! Pack send buffer, send east slice to update west boundary of 
  ! the east neighbor
  !
  DO i=1,bdrwidth
    tem1(i) = xs(nx-3-bdrwidth+i)
  END DO

  CALL mpi_sendrecv(tem1,bdrwidth,MPI_REAL,dest,  mptag+tag_w,       &
                    tem2,bdrwidth,MPI_REAL,source,mptag+tag_w,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update West boundary data
  !
  IF ( loc_x /= 1 ) THEN
    DO i=1,bdrwidth
      xsext(i-bdrwidth) = tem2(i)
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  ! 
  ! receive from
  !
  IF(loc_x == nproc_x) THEN        ! Last processor in a row
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

!---------------------------- U-GRID -----------------------------------

  !
  ! Pack send buffer, west slice for east boundary of 
  ! the west neighbor
  !
  DO i=1,bdrwidth
    tem1(i) = x(3+i)
  END DO

  CALL mpi_sendrecv(tem1,bdrwidth,MPI_REAL,dest,  mptag+tag_e,       &
                    tem2,bdrwidth,MPI_REAL,source,mptag+tag_e,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x )  THEN
    DO i=1,bdrwidth
      xext(nx+i) = tem2(i)
    END DO
  END IF

!---------------------------- MASS -----------------------------------

  !
  ! Pack send buffer, west slice for east boundary of 
  ! the west neighbor
  !
  DO i=0,bdrwidth
    tem1(i+1) = xs(3+i)
  END DO

  CALL mpi_sendrecv(tem1,bdrwidth+1,MPI_REAL,dest,  mptag+tag_e,       &
                    tem2,bdrwidth+1,MPI_REAL,source,mptag+tag_e,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x )  THEN
    DO i=0,bdrwidth
      xsext(nx+i) = tem2(i+1)
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == 1) THEN             ! the south most processor in a column
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  ! 
  ! receive from
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor in a column
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+nproc_x*loc_y)
  END IF

!---------------------------- V-GRID -----------------------------------

  !
  ! Pack send buffer, send south slice to update north boundary of 
  ! the south neighbor
  !
  DO j=1,bdrwidth
    tem1(j) = y(3+j)
  END DO

  CALL mpi_sendrecv(tem1,bdrwidth,MPI_REAL,dest,  mptag+tag_n,       &
                    tem2,bdrwidth,MPI_REAL,source,mptag+tag_n,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y )  THEN
    DO j=1,bdrwidth
      yext(ny+j) = tem2(j)
    END DO
  END IF

!---------------------------- MASS -----------------------------------

  !
  ! Pack send buffer, send south slice to update north boundary of 
  ! the south neighbor
  !
  DO j=0,bdrwidth
    tem1(j+1) = ys(3+j)
  END DO

  CALL mpi_sendrecv(tem1,bdrwidth+1,MPI_REAL,dest,  mptag+tag_n,       &
                    tem2,bdrwidth+1,MPI_REAL,source,mptag+tag_n,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y )  THEN
    DO j=0,bdrwidth      
      ysext(ny+j) = tem2(j+1)
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+nproc_x*loc_y)
  END IF

  ! 
  ! receive from
  !
  IF(loc_y == 1) THEN            ! The south most processor
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+nproc_x*(loc_y-2))
  END IF

!---------------------------- V-GRID -----------------------------------

  !
  ! Pack send buffer, north slice for south boundary of 
  ! the north neighbor
  !
  DO j=1,bdrwidth
    tem1(j) = y(ny-3-bdrwidth+j)
  END DO

  CALL mpi_sendrecv(tem1,bdrwidth,MPI_REAL,dest,  mptag+tag_s,       &
                    tem2,bdrwidth,MPI_REAL,source,mptag+tag_s,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update south boundary data
  !
  IF ( loc_y /= 1 )  THEN
    DO j=1,bdrwidth
      yext(j-bdrwidth) = tem2(j)
    END DO
  END IF


!---------------------------- MASS -----------------------------------

  !
  ! Pack send buffer, north slice for south boundary of 
  ! the north neighbor
  !
  DO j=1,bdrwidth
    tem1(j) = ys(ny-3-bdrwidth+j)
  END DO

  CALL mpi_sendrecv(tem1,bdrwidth,MPI_REAL,dest,  mptag+tag_s,       &
                    tem2,bdrwidth,MPI_REAL,source,mptag+tag_s,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update south boundary data
  !
  IF ( loc_y /= 1 )  THEN
    DO j=1,bdrwidth
      ysext(j-bdrwidth) = tem2(j)
    END DO
  END IF


!---------------------------- CLOSING SUBROUTINE -----------------------

  DEALLOCATE(tem1,tem2)
#endif

  RETURN
END SUBROUTINE set_coord_halo
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE EXCHANGE_HALO               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE exchange_halo(svarin,nx,ny,nz,                               &
                         varout,mxs,mxe,mys,mye,istagger,istatus)

!#######################################################################
!  Expand memory array for message passing purpose

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx, ny, nz
  REAL,    INTENT(IN) :: svarin(nx,ny,nz)
  INTEGER, INTENT(IN) :: mxs, mxe, mys, mye
  INTEGER, INTENT(IN) :: istagger
  REAL,    INTENT(INOUT) :: varout(mxs:mxe,mys:mye,nz)
  INTEGER, INTENT(OUT)   :: istatus
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

  REAL, PARAMETER :: MISSING_DATA = HUGE(1.0)

#ifdef MPI_PARALLEL
  INCLUDE 'mpif.h'
  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: imstat
#endif

  INTEGER :: source, dest
  INTEGER :: mptag           ! Unique MPI id used for this BC update.
                             ! mptag + tag_w for west boundary
                             ! mptag + tag_e for east boundary

  INTEGER :: i, j, k, n
  INTEGER :: si, sj

  INTEGER :: bdrwidth, msize

  REAL, ALLOCATABLE :: tem1(:,:,:), tem2(:,:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  bdrwidth = 1-mxs

  IF (bdrwidth <= 0) RETURN

  si = 0
  sj = 0
  IF (istagger == 1) si = 1
  IF (istagger == 2) sj = 1

#ifdef MPI_PARALLEL

  IF (loc_x == 1) THEN
    DO k = 1,nz
      DO j = mys,mye
        DO i = mxs,0
          varout(i,j,k) = MISSING_DATA
        END DO
      END DO
    END DO
  END IF

  IF (loc_x == nproc_x) THEN
    DO k = 1,nz
      DO j = mys,mye
        DO i = nx+si,mxe
          varout(i,j,k) = MISSING_DATA
        END DO
      END DO
    END DO
  END IF

  IF (loc_y == 1) THEN
    DO k = 1,nz
      DO j = mys,0
        DO i = mxs,mxe
          varout(i,j,k) = MISSING_DATA
        END DO
      END DO
    END DO
  END IF

  IF (loc_y == nproc_y) THEN
    DO k = 1,nz
      DO j = ny+sj,mye
        DO i = mxs,mxe
          varout(i,j,k) = MISSING_DATA
        END DO
      END DO
    END DO
  END IF

  CALL inctag
  mptag = gentag

  ALLOCATE(tem1(bdrwidth,ny,nz), STAT = istatus)
  ALLOCATE(tem2(bdrwidth,ny,nz), STAT = istatus)

!-----------------------------------------------------------------------
!
!  Set the west boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == nproc_x) THEN       ! last processor in a row
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  ! 
  ! receive from
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, send east slice to update west boundary of 
  ! the east neighbor
  !
  DO k = 1,nz
    DO j = 1,ny
      DO i = 1, bdrwidth
          n = nx-3-bdrwidth+i
          tem1(i,j,k) = svarin(n,j,k)
      END DO
    END DO
  END DO

  CALL mpi_sendrecv(tem1,bdrwidth*ny*nz,MPI_REAL,dest,  mptag+tag_w,       &
                    tem2,bdrwidth*ny*nz,MPI_REAL,source,mptag+tag_w,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update West boundary data
  !
  IF ( loc_x /= 1 ) THEN
    DO k = 1,nz
      DO j = 1,ny
        DO i=1,bdrwidth
           varout(i-bdrwidth,j,k) = tem2(i,j,k)
        END DO
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Set the east boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_x == 1) THEN             ! First processor in a row
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x-1+nproc_x*(loc_y-1))
  END IF

  ! 
  ! receive from
  !
  IF(loc_x == nproc_x) THEN        ! Last processor in a row
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+1+nproc_x*(loc_y-1))
  END IF

  !
  ! Pack send buffer, west slice for east boundary of 
  ! the west neighbor
  !
  DO k = 1,nz
    DO j = 1,ny
      DO i = 1,bdrwidth
        tem1(i,j,k) = svarin(2+i+si,j,k)
      END DO
    END DO
  END DO

  CALL mpi_sendrecv(tem1,bdrwidth*ny*nz,MPI_REAL,dest,  mptag+tag_e,       &
                    tem2,bdrwidth*ny*nz,MPI_REAL,source,mptag+tag_e,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update east boundary data
  !
  IF ( loc_x /= nproc_x )  THEN
    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,bdrwidth
          varout(nx-1+si+i,j,k) = tem2(i,j,k)
        END DO
      END DO
    END DO
  END IF


  DEALLOCATE(tem1,tem2)

  ALLOCATE(tem1(mxs:mxe,bdrwidth,nz), STAT = istatus)
  ALLOCATE(tem2(mxs:mxe,bdrwidth,nz), STAT = istatus)

  msize = (mxe-mxs+1)*bdrwidth*nz

!-----------------------------------------------------------------------
!
!  Set the north boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == 1) THEN             ! the south most processor in a column
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  ! 
  ! receive from
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor in a column
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+nproc_x*loc_y)
  END IF

  !
  ! Pack send buffer, send south slice to update north boundary of 
  ! the south neighbor
  !
  DO k = 1,nz
    DO j = 1,bdrwidth
      DO i = mxs,mxe
        tem1(i,j,k) = varout(i,2+j+sj,k)
      END DO
    END DO
  END DO

  CALL mpi_sendrecv(tem1,msize,MPI_REAL,dest,  mptag+tag_n,       &
                    tem2,msize,MPI_REAL,source,mptag+tag_n,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update north boundary data
  !
  IF ( loc_y /= nproc_y )  THEN
    DO k = 1,nz
      DO j = 1,bdrwidth
        DO i = mxs,mxe
          varout(i,ny-1+j+sj,k) = tem2(i,j,k)
        END DO
      END DO
    END Do
  END IF

!-----------------------------------------------------------------------
!
!  Set the south boundary conditions
!
!-----------------------------------------------------------------------

  !
  ! send destination
  !
  IF(loc_y == nproc_y) THEN       ! The north most processor
    dest = MPI_PROC_NULL
  ELSE
    dest = proc(loc_x+nproc_x*loc_y)
  END IF

  ! 
  ! receive from
  !
  IF(loc_y == 1) THEN            ! The south most processor
    source = MPI_PROC_NULL
  ELSE
    source = proc(loc_x+nproc_x*(loc_y-2))
  END IF

  !
  ! Pack send buffer, north slice for south boundary of 
  ! the north neighbor
  !
  DO k = 1, nz
    DO j = 1,bdrwidth
      DO i = mxs,mxe
        tem1(i,j,k) = varout(i,ny-3-bdrwidth+j,k)
      END DO
    END DO
  END DO

  CALL mpi_sendrecv(tem1,msize,MPI_REAL,dest,  mptag+tag_s,       &
                    tem2,msize,MPI_REAL,source,mptag+tag_s,       &
                    MPI_COMM_WORLD,mpi_status,imstat)

  !
  ! Unpack receive buffer, update south boundary data
  !
  IF ( loc_y /= 1 )  THEN
    DO k = 1,nz
      DO j = 1,bdrwidth
        DO i = mxs,mxe
          varout(i,j-bdrwidth,k) = tem2(i,j,k)
        END DO
      END Do
    END DO
  END IF

#endif

!---------------------------- CLOSING SUBROUTINE -----------------------

  RETURN
END SUBROUTINE exchange_halo
