
SUBROUTINE pltall(ndim)
!*********************************************************************
! plot ndim-variables in x-y planes for all grids
!*********************************************************************

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'

!*********************************************************************
!  plot x-y planes of ndim-D variables on grid mptr
!*********************************************************************

  CALL xcltyp(0)
  CALL xcfull
  CALL xhilit(0)
  mptr=1

!*********************************************************************
!  goto 44
! plot vector component ivar=2 u
!*********************************************************************

  ivar = 2
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL pltvxy(ndim,2,5,11, 1, 2)
  CALL xframe
  CALL pltvxy(ndim,2,5,11, 1, nzz/2)
  CALL xframe

!*********************************************************************
! plot vector component ivar=5 v
!*********************************************************************

  ivar = 5
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL pltvxy(ndim,2,5,11, 2, 2)
  CALL xframe
  CALL pltvxy(ndim,2,5,11, 2, nzz/2)
  CALL xframe

!*********************************************************************
! plot scalar w
!*********************************************************************

  ivar=8
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL pltsxy(ndim,ivar,2)
  CALL xframe
  CALL pltsxy(ndim,ivar,nzz/2)
  CALL xframe

!*********************************************************************
! plot scalar ptprt
!*********************************************************************

  ivar=11
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL pltsxy(ndim,ivar,2)
  CALL xframe
  CALL pltsxy(ndim,ivar,nzz/2)
  CALL xframe

!*********************************************************************
! plot scalar pprt
!*********************************************************************

  ivar=14
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL pltsxy(ndim,ivar,2)
  CALL xframe
  CALL pltsxy(ndim,ivar,nzz/2)
  CALL xframe

!*********************************************************************
! plot scalar qvprt
!*********************************************************************

  ivar=17
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL pltsxy(ndim,ivar,2)
  CALL xframe
  CALL pltsxy(ndim,ivar,nzz/2)
  CALL xframe

!*********************************************************************
! plot scalar qcprt
!*********************************************************************

  ivar=20
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL pltsxy(ndim,ivar,2)
  CALL xframe
  CALL pltsxy(ndim,ivar,nzz/2)
  CALL xframe

!*********************************************************************
! plot scalar qrprt
!*********************************************************************

  ivar=23
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL pltsxy(ndim,ivar,2)
  CALL xframe
  CALL pltsxy(ndim,ivar,nzz/2)
  CALL xframe

!*********************************************************************

  RETURN
END SUBROUTINE pltall

!*********************************************************************
!*********************************************************************
!*********************************************************************



SUBROUTINE gtrcrd(mptr,ndim,ivar,x,y,nx1,ny1)
!*********************************************************************
! To return the coordinate arrays x and y of variable ivar on grid mptr
! relative to its own grid origin
! nx1, ny1 give the dimension of the portion of array actually used.
!*********************************************************************

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  REAL :: x(1), y(1)

  CALL chkdim(ndim,'GTRCRD')

  nx=node(5,mptr)
  ny=node(6,mptr)
  cosr=rnode(21,mptr)
  sinr=rnode(22,mptr)
  dx =rnode(9 ,mptr)
  dy =rnode(10,mptr)
  xor=rnode(1,mptr)
  yor=rnode(2,mptr)

  IF(ndim == 2) THEN
    IF(ivar > nxy2d)THEN
      WRITE(6,5)  ivar,nxy2d
      WRITE(6,'(1x,a)') 'Job stopped in GTRCRD!'
      STOP
    END IF
    xshift=stgxy(1,ivar)
    yshift=stgxy(2,ivar)
    nx1=nx-idmxy(1,ivar)
    ny1=ny-idmxy(2,ivar)
  ELSE
    xshift=stgxyz(1,ivar)
    yshift=stgxyz(2,ivar)
    nx1=nx-idmxyz(1,ivar)
    ny1=ny-idmxyz(2,ivar)
    IF(ivar > nxyz3d)THEN
      WRITE(6,10) ivar,nxyz3d
      WRITE(6,'(1x,a)') 'Job stopped in GTRCRD!'
      STOP
    END IF
  END IF

  DO j=1,ny
    DO i=1,nx
      ij=i+(j-1)*nx
      x(ij)= (i-1+xshift)*dx
      y(ij)= (j-1+yshift)*dy
    END DO
  END DO
  5     FORMAT                                                          &
         ('Error in GTACRD: variable number for 2-d x-y array ivar='    &
              ,i3,' > nxy2d=',i3)
  10    FORMAT('Error in GTACRD: variable number for 3-d array ivar=',i3, &
                 ' > nxyz3d=',i3)
  RETURN
END SUBROUTINE gtrcrd

!*********************************************************************
!*********************************************************************
!*********************************************************************


SUBROUTINE gtacrd(mptr,ndim,ivar,x,y,nx1,ny1)
!*********************************************************************
! To return the absolute coordinate arrays x and y of variable ivar
! on grid mptr relative to the absolute origin.
! nx1, ny1 give the dimension of the portion of array actually used.
!*********************************************************************

  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  REAL :: x(1),y(1)

  CALL chkdim(ndim,'GTRCRD')

  nx=node(5,mptr)
  ny=node(6,mptr)
  cosr=rnode(21,mptr)
  sinr=rnode(22,mptr)
  dx =rnode(9 ,mptr)
  dy =rnode(10,mptr)
  xor=rnode(1,mptr)
  yor=rnode(2,mptr)

  IF(ndim == 2) THEN
    IF(ivar > nxy2d)THEN
      WRITE(6,5)  ivar,nxy2d
      WRITE(6,'(1x,a)') 'Job stopped in GTACRD!'
      STOP
    END IF
    xshift=stgxy(1,ivar)
    yshift=stgxy(2,ivar)
    nx1=nx-idmxy(1,ivar)
    ny1=ny-idmxy(2,ivar)
  ELSE
    IF(ivar > nxyz3d)THEN
      WRITE(6,10) ivar,nxyz3d
      WRITE(6,'(1x,a)') 'Job stopped in GTACRD!'
      STOP
    END IF
    xshift=stgxyz(1,ivar)
    yshift=stgxyz(2,ivar)
    nx1=nx-idmxyz(1,ivar)
    ny1=ny-idmxyz(2,ivar)
  END IF

  DO j=1,ny
    DO i=1,nx
      xx = (i-1+xshift)*dx
      yy = (j-1+yshift)*dy
      ij=i+(j-1)*nx
      x(ij) = xor+xx*cosr-yy*sinr
      y(ij) = yor+xx*sinr+yy*cosr
    END DO
  END DO
  5     FORMAT                                                          &
           ('Error in GTACRD: variable number for 2-d x-y array ivar=', &
           i3,' > nxy2d=',i3)
  10    FORMAT('Error in GTACRD: variable number for 3-d array ivar=',  &
              i3,' > nxyz3d=',i3)
  RETURN
END SUBROUTINE gtacrd

!*********************************************************************
!*********************************************************************
!*********************************************************************



SUBROUTINE gtxypl( a,tmp,m,n,l,k )
!*********************************************************************
! copy a x-y plane from a 3-d field
!*********************************************************************

  DIMENSION a(m,n,l),tmp(m,n)
  DO j=1,n
    DO i=1,m
      tmp(i,j) = a(i,j,k)
    END DO
  END DO
  RETURN
END SUBROUTINE gtxypl

!*********************************************************************
!*********************************************************************
!*********************************************************************


  INTEGER FUNCTION nxtgrd(mptr0)
!*********************************************************************
! return a grid pointer that is next to mptr0 on the same level or the
! first one on the level that is one-level finer than the level of mptr0
!*********************************************************************


  INCLUDE 'nodal.inc'

  mptr=mptr0
  IF(mptr == 0) RETURN
  level=node(4,mptr)
  5     mptr=node(10,mptr)
  IF(mptr /= 0) GO TO 20
  level = level+1
  IF(level > lfine)THEN
    nxtgrd=0
    GO TO 999
  END IF
  nxtgrd=lstart(level)
  GO TO 999
  20    nxtgrd=mptr
  999   RETURN
  END FUNCTION nxtgrd

!*********************************************************************
!*********************************************************************
!*********************************************************************


  INTEGER FUNCTION lstgrd(mptr0)
!*********************************************************************
! return a grid pointer that is next to mptr0 on the same level or the
! first one on the level that is one-level coarser than the level of mptr0
!*********************************************************************

  INCLUDE 'nodal.inc'

  mptr=mptr0
  IF(mptr == 0) RETURN
  level=node(4,mptr)
  5     mptr=node(10,mptr)
  IF(mptr /= 0) GO TO 20
  level = level-1
  IF(level <= 0)THEN
    lstgrd=0
    GO TO 999
  END IF
  lstgrd=lstart(level)
  GO TO 999
  20    lstgrd=mptr
  999   RETURN
  END FUNCTION lstgrd

!*********************************************************************
!*********************************************************************
!*********************************************************************


SUBROUTINE pltsxy(ndim,ivar, kk)
!*********************************************************************
! plot x-y slices of ndim-D scalar ivar at k=kk.
!*********************************************************************
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  CHARACTER (LEN=132) :: title

  pi=4.0*ATAN(1.0)
  CALL xqmap (x1,x2,y1,y2)
  CALL xqlmsk(nmask)

!*********************************************************************
! plot a over all frame and axes.
!*********************************************************************

  mptr=mstart
  nx=node(5,mptr)
  ny=node(6,mptr)
  nz=node(14,mptr)
  pi=4.0*ATAN(1.0)
  xor=rnode(1, mptr)
  yor=rnode(2, mptr)
  dx=rnode(9 , mptr)
  dy=rnode(10, mptr)
  xl=(nx-1)*dx
  yl=(ny-1)*dy
  IF(yl > xl) THEN
    xm = (1.0-0.8*xl/yl)/2
    CALL xpspac( xm, 1.0-xm, 0.1, 0.9)
  ELSE
    ym = (1.0-0.8*yl/xl)/2
    CALL xpspac(0.1,0.9, ym, 1.0-ym)
  END IF

  CALL xmap  (0.0, xl, 0.0, yl)
  CALL xqpspc(p1,p2,p3,p4)
  CALL xqmap(xa,xb,ya,yb)
  CALL xaxsca(0.0, xl, dx, 0.0, yl, dy)
  WRITE(title,'(i1,''-D VAR NO.'',I2,'' K='',I3)')ndim,ivar,kk
  CALL xchsiz(0.025*yl)
  CALL xcharc(0.5*xl, 1.05*yl, title(1:21))
  PRINT*
  PRINT*,'Plotting ',title(1:21), ' for grid ', mptr

!*********************************************************************
! loop over all grids
!*********************************************************************

  mptr = lstart(lfine)
  100   CONTINUE
  nx=node(5,mptr)
  ny=node(6,mptr)
  nz=node(14,mptr)

  CALL xqlmsk(level)
  CALL xunmsk(1)
  IF(mptr /= mstart) THEN
    CALL xpenup(rnode(7,mptr),rnode(8,mptr))
    DO i=1,7,2
      CALL xpendn(rnode(i,mptr),rnode(i+1,mptr))
    END DO
  END IF
  CALL xrsmsk(level)
  IF(ndim == 3) in1 = igtxyz(mptr,ivar,1 )
  IF(ndim == 2) in1 = igtnxy(mptr,ivar,1 )
  n3rd=nz
  IF(ndim == 2) n3rd=1
!
  needsp = nx*ny
  ntem1=igetsp(needsp)
  ntem2=igetsp(needsp)
  ntem3=igetsp(needsp)
  CALL gtxypl(a(in1),a(ntem1),nx,ny,n3rd,MAX(1,kk))
  CALL gtacrd(mptr,ndim,ivar,a(ntem2),a(ntem3),nxa,nya)
  IF(mptr == lstart(lfine) ) THEN
    zinc=1.0
    mode=1
  ELSE
    mode=2  ! now use the zinc of first grid lstart(lfine).
  END IF

  IF(mptr == mstart)THEN
    CALL xconts2(a(ntem1),a(ntem2),a(ntem3),nx,nxa,nya,.true.,          &
         zinc,mode)
  ELSE
    CALL xconts2(a(ntem1),a(ntem2),a(ntem3),nx,nxa,nya,.false.,         &
         zinc,mode)
    xor=rnode(1, mptr)
    yor=rnode(2, mptr)
    dx =rnode(9 , mptr)
    dy =rnode(10, mptr)
    xl=(nx-1)*dx
    yl=(ny-1)*dy
    cosr=rnode(21,mptr)
    sinr=rnode(22,mptr)
    ang = ASIN (sinr )*180/pi
    xo = xor+0.5*(dx*cosr-dy*sinr)
    yo = yor+0.5*(dx*sinr+dy*cosr)
    CALL xrmask(xo,yo,xl-dx,yl-dx,ang)
  END IF

  CALL resett
  mptr= lstgrd(mptr)
  IF(mptr /= 0)GO TO 100

  CALL xmap(x1,x2,y1,y2)
  CALL xrsmsk(nmask)
  RETURN
END SUBROUTINE pltsxy

!*********************************************************************
!*********************************************************************
!*********************************************************************


SUBROUTINE pltvxy(ndim,ivar1,ivar2,iscalr, kk, iplot)
!*********************************************************************
! plot contour map of one of two components (ivar1, ivar2) of
! a ndim-D vector projected to absolute coordinate system, i.e.
! that of the base grid. The coordinates for scalar variable iscalr
! are used.
!
! iplot =1, ivar1 isplotted, else ivar2 is plotted.
!*********************************************************************
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  CHARACTER (LEN=132) :: title

  pi=4.0*ATAN(1.0)
  CALL xqmap (x1,x2,y1,y2)
  CALL xqlmsk(nmask)

  IF(iplot == 1) THEN
    ivar=ivar1
  ELSE
    ivar=ivar2
  END IF

  mptr=mstart
  nx=node(5,mptr)
  ny=node(6,mptr)
  nz=node(14,mptr)
  pi=4.0*ATAN(1.0)
  xor=rnode(1, mptr)
  yor=rnode(2, mptr)
  dx=rnode(9 , mptr)
  dy=rnode(10, mptr)
  xl=(nx-1)*dx
  yl=(ny-1)*dy
  IF(yl > xl) THEN
    xm = (1.0-0.8*xl/yl)/2
    CALL xpspac( xm, 1.0-xm, 0.1, 0.9)
  ELSE
    ym = (1.0-0.8*yl/xl)/2
    CALL xpspac(0.1,0.9, ym, 1.0-ym)
  END IF

  CALL xmap  (0.0, xl, 0.0, yl)
  CALL xqpspc(p1,p2,p3,p4)
  CALL xqmap(xa,xb,ya,yb)
  CALL xaxsca(0.0, xl, dx, 0.0, yl, dy)

  WRITE(title,'(i1,''-D VAR NO.'',I2,'' K='',I3)')ndim,ivar,kk
  CALL xchsiz(0.025*yl)
  CALL xcharc(0.5*xl, 1.05*yl, title(1:21))
  PRINT*
  PRINT*,'Plotting ',title(1:21), ' for grid ', mptr

  mptr = lstart(lfine)
  100   CONTINUE
  nx=node(5,mptr)
  ny=node(6,mptr)
  nz=node(14,mptr)

  CALL xqlmsk(level)
  CALL xunmsk(1)
  IF(mptr /= mstart) THEN
    CALL xpenup(rnode(7,mptr),rnode(8,mptr))
    DO i=1,7,2
      CALL xpendn(rnode(i,mptr),rnode(i+1,mptr))
    END DO
  END IF
  CALL xrsmsk(level)
  IF(ndim == 3) in1 = igtxyz(mptr,ivar1,1 )
  IF(ndim == 2) in1 = igtnxy(mptr,ivar1,1 )
  IF(ndim == 3) in2 = igtxyz(mptr,ivar2,1 )
  IF(ndim == 2) in2 = igtnxy(mptr,ivar2,1 )
  n3rd=nz
  IF(ndim == 2) n3rd=1
!
  needsp = nx*ny
  ix=igetsp(needsp)
  iy=igetsp(needsp)
  iu=igetsp(needsp)
  iv=igetsp(needsp)
  CALL tranuv(mptr,a(in1),a(in2),a(iu),a(iv),nx,ny,n3rd,MAX(1,kk))
  CALL gtacrd(mptr,ndim,iscalr,a(ix),a(iy),nxa,nya)
  IF(mptr == lstart(lfine) ) THEN
    zinc=1.0
    mode=1
  ELSE
    mode=2  ! now use the zinc of first grid lstart(lfine).
  END IF

  IF(iplot == 1) THEN
    iptr=iu
  ELSE
    iptr=iv
  END IF

  IF(mptr == mstart)THEN
    CALL xconts2(a(iptr),a(ix),a(iy),nx,nxa,nya,.true.,                 &
         zinc,mode)
  ELSE
    CALL xconts2(a(iptr),a(ix),a(iy),nx,nxa,nya,.false.,                &
         zinc,mode)
    xor=rnode(1, mptr)
    yor=rnode(2, mptr)
    dx =rnode(9 , mptr)
    dy =rnode(10, mptr)
    xl=(nx-1)*dx
    yl=(ny-1)*dy
    cosr=rnode(21,mptr)
    sinr=rnode(22,mptr)
    ang = ASIN (sinr )*180/pi
    xo = xor+0.5*(dx*cosr-dy*sinr)
    yo = yor+0.5*(dx*sinr+dy*cosr)
    CALL xrmask(xo,yo,xl-dx,yl-dx,ang)
  END IF

  CALL resett
  mptr= lstgrd(mptr)

  IF(mptr /= 0)GO TO 100
  CALL xmap(x1,x2,y1,y2)
  CALL xrsmsk(nmask)
  RETURN
END SUBROUTINE pltvxy



!*********************************************************************
!*********************************************************************
!*********************************************************************





SUBROUTINE pltgrd(mptr, ndim)
!*********************************************************************
! plot contour maps of ndim-D fields on grid mptr. Only selectedx-y slices
! are plotted.
!*********************************************************************

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'

!*********************************************************************
!  plot x-y planes of ndim-D variables on grid mptr
!
!    goto 111
!*********************************************************************

  ivar=1
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL plotxy(mptr,ndim,ivar,1)
  CALL xframe
  CALL plotxy(mptr,ndim,ivar,nzz)
  CALL xframe
!
  ivar=3
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL plotxy(mptr,ndim,ivar,1)
  CALL xframe
  CALL plotxy(mptr,ndim,ivar,nzz)
  CALL xframe

  ivar=5
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL plotxy(mptr,ndim,ivar,1)
  CALL xframe
  CALL plotxy(mptr,ndim,ivar,nzz)
  CALL xframe

  111     CONTINUE
  ivar=7
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL plotxy(mptr,ndim,ivar,1)
  CALL xframe
  CALL plotxy(mptr,ndim,ivar,nzz)
  CALL xframe

  RETURN
  ivar=9
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL plotxy(mptr,ndim,ivar,1)
  CALL xframe
  CALL plotxy(mptr,ndim,ivar,nzz)
  CALL xframe

  RETURN
  ivar=11
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL plotxy(mptr,ndim,ivar,1)
  CALL xframe
  CALL plotxy(mptr,ndim,ivar,nzz)
  CALL xframe

  ivar=12
  nzz = node(14,mptr)-idmxyz(3,ivar)
  CALL plotxy(mptr,ndim,ivar,1)
  CALL xframe
  CALL plotxy(mptr,ndim,ivar,nzz)
  CALL xframe
  RETURN
END SUBROUTINE pltgrd

!*********************************************************************
!*********************************************************************
!*********************************************************************


SUBROUTINE plotxy(mptr,ndim,ivar, kk)
!*********************************************************************
! plot x-y slices ofr ndim-d field ivar on grid mptr.
! no coordinate transformation is done.
!*********************************************************************

  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  CHARACTER (LEN=132) :: title

!
  nx=node(5,mptr)
  ny=node(6,mptr)
  nz=node(14,mptr)

  pi=4.0*ATAN(1.0)
  xor=rnode(1, mptr)
  yor=rnode(2, mptr)
  dx=rnode(9 , mptr)
  dy=rnode(10, mptr)
  xl=(nx-1)*dx
  yl=(ny-1)*dy

  CALL xqmap (x1,x2,y1,y2)

  IF(yl > xl) THEN
    CALL xpspac(0.1+(yl-xl)/(yl*2),0.9-(yl-xl)/(yl*2),0.1,0.8)
  ELSE
    CALL xpspac(0.1,0.9,0.1+(xl-yl)/(xl*2),0.9-(xl-yl)/(xl*2))
  END IF
  CALL xmap  (0.0, xl, 0.0, yl)
  CALL xaxsca(0.0, xl, dx, 0.0, yl, dy)

  IF(ndim == 3) in1 = igtxyz(mptr,ivar,1 )
  IF(ndim == 2) in1 = igtnxy(mptr,ivar,1 )
  n3rd=nz
  IF(ndim == 2) n3rd=1
!
  needsp = nx*ny
  ntem1=igetsp(needsp)
  ntem2=igetsp(needsp)
  ntem3=igetsp(needsp)
  CALL gtxypl(a(in1),a(ntem1),nx,ny,n3rd,MAX(1,kk))
!  call gtacrd(mptr,ndim,ivar,a(ntem2),a(ntem3),nxa,nya)
  CALL gtrcrd(mptr,ndim,ivar,a(ntem2),a(ntem3),nxa,nya)
  CALL xconts1(a(ntem1),a(ntem2),a(ntem3),nx,nxa,nya,.true.)
  WRITE(title,'(''VAR NO.'',I2,'' K='',I3,'' FOR GRID '',I2)')          &
        ivar, kk, mptr
  CALL xchsiz(0.025*yl)
  CALL xcharc(0.5*xl, 1.05*yl, title(1:27))

  CALL xmap  (x1,x2,y1,y2)
  WRITE(title,'(''X-Y PLANE FOR K='',I2,'' OF '',                       &
      i1,''-d var. no.'',i3,'' of grid'',i2)') kk,ndim,ivar,mptr
  PRINT*,title(1:60)
!  call WRIGAR(A(ntem1),1,nx,1,ny,1,1,1,nx,1,ny,1,1,
!    :     title(1:60),0.0 ,1)
  CALL resett
  RETURN
END SUBROUTINE plotxy


!*********************************************************************
!*********************************************************************
!*********************************************************************

SUBROUTINE xconts1(z,x,y,md,m,n, limit)
!*********************************************************************
! plot a contour map for field Z.
! if limit=.false., max, min, inc are not plotted.
!*********************************************************************

  REAL :: z(md,n),x(md,n),y(md,n),cl(150)
  INTEGER :: iwrk(10000)
  LOGICAL :: limit

  IF(m*n > 10000)THEN
    PRINT*,'Working array IW in plotting routine xconts1',              &
        ' defined too small'
    WRITE(6,'(a,i6,a,i6)')'Space needed:',m*n,', allocated:',10000
    PRINT*,'No further plotting performed in this routine.'
    RETURN
  END IF

  zmin=z(1,1)
  zmax=zmin
  DO j=1,n
    DO i=1,m
      zmax=MAX(zmax,z(i,j))
      zmin=MIN(zmin,z(i,j))
    END DO
  END DO
  PRINT*,'zmax, zmin ', zmax, zmin
  cl(1)=0.0
  IF(zmax-zmin > 1.0E-8)THEN
    cl(2)=cl(1)+ xfinc(zmax-zmin)
    IF(cl(2)-cl(1) == 0.0) cl(2)=cl(1)+1.0
    mode=1
    CALL xconta(z,x,y,iwrk,md,m,n,cl,ncl,mode)
  ELSE
    cl(2)=1.0
    ncl=2
  END IF
  IF(.NOT.limit) RETURN
!  zmax =cl(ncl)
!  zmin =cl( 1)
  zinc=cl(MIN(2,ncl))-cl(1)
  CALL xclimt1(zmax,zmin,zinc,0.0)
  RETURN
END SUBROUTINE xconts1


!*********************************************************************
!*********************************************************************
!*********************************************************************

SUBROUTINE xconts2(z,x,y,md,m,n, limit, zinc0, model)
!*********************************************************************
! produce contour plot for field Z with coordinate X(i,j),Y(i,j).
! limit = .t. then zmin,zmax & zinc are plotted above the top
! model = 1, zinc is determined inside and on output zinc0 = zinc.
!    = 2, zinc = zinc0.
!*********************************************************************

  REAL :: z(md,n),x(md,n),y(md,n),cl(150)
  INTEGER :: iwrk(10000)
  LOGICAL :: limit

  IF(m*n > 10000)THEN
    PRINT*,'Working array IW in plotting routine xconts1',              &
        ' defined too small'
    WRITE(6,'(a,i6,a,i6)')'Space needed:',m*n,', allocated:',10000
    PRINT*,'No further plotting performed in this routine.'
    RETURN
  END IF

  zmin=z(1,1)
  zmax=zmin
  DO j=1,n
    DO i=1,m
      zmax=MAX(zmax,z(i,j))
      zmin=MIN(zmin,z(i,j))
    END DO
  END DO
  PRINT*,'zmax, zmin ', zmax, zmin
  cl(1)=0.0
  IF(zmax-zmin > 1.0E-8)THEN
    IF(model /= 2)THEN
      cl(2)=cl(1)+ xfinc(zmax-zmin)
      IF(cl(2)-cl(1) == 0.0) cl(2)=cl(1)+1.0
      mode=1
    ELSE
      cl(2)=zinc0
      mode=2
    END IF
    CALL xconta(z,x,y,iwrk,md,m,n,cl,ncl,mode)
  ELSE
    cl(2)=1.0
    ncl=2
  END IF
!  zmax =cl(ncl)
!  zmin =cl( 1)
  zinc=cl(MIN(2,ncl))-cl(1)
  zinc0 = zinc
  IF(.NOT.limit) RETURN
  CALL xclimt1(zmax,zmin,zinc,0.0)
  RETURN
END SUBROUTINE xconts2


!*********************************************************************
!*********************************************************************
!*********************************************************************

  REAL FUNCTION xfinc(x)
!*********************************************************************
! to automatically divide domain (0,x) to a number of subdomain
! with intervale xfinc which is >=4 and =<16 for fold=1.0)
!*********************************************************************

  ipower= INT(  LOG10(x) )
  d= INT(x/(10.0**ipower))
  fold=1.0
  xfinc=0.1*x
  IF(d >= 0.0.AND.d < 3.0) THEN
    xfinc=2.0*10.0**(ipower-1)
  ELSE IF(d >= 3.0.AND.d < 7.0) THEN
    xfinc=5.0*10.0**(ipower-1)*fold
  ELSE IF(d >= 7.0.AND.d < 10.) THEN
    xfinc=1.0*10.0** ipower*fold
  END IF
  IF(xfinc == 0.0) xfinc=x*0.1
  RETURN
  END FUNCTION xfinc

!*********************************************************************
!*********************************************************************
!*********************************************************************

SUBROUTINE xclimt1(fmax,fmin,finc )
!
  CHARACTER (LEN=150) :: ch
  CALL xqmap(xl,xr,yb,yt)
  CALL xqrang( xrg, yrg )
  CALL xqchmg( siz0 )
  siz=0.02*SQRT( xrg*yrg)
  CALL xchmag(siz)
  WRITE(ch,'(''MIN='',G9.3E2,'' MAX='',G9.3E2,                          &
      '' inc='',g9.3E2)')fmin,fmax,finc
  CALL xcharc( 0.5*(xl+xr), yt+0.02*(yt-yb),ch(1:41) )
  CALL xchmag( siz0)
  RETURN
END SUBROUTINE xclimt1

!*********************************************************************
!*********************************************************************
!*********************************************************************

SUBROUTINE tranuv(mptr,u,v,uc,vc,nx,ny,nz,kk)
!*********************************************************************
! transform vector field (u,v) on grid mptr to absolute coordinate
! and copy x-y slice at k=kk to 2-D field (uc,vc).
! these vectors are located in the center of staggered grid cells.
!*********************************************************************

  REAL :: u(nx,ny,nz),v(nx,ny,nz), uc(nx,ny),vc(nx,ny)
  INCLUDE 'nodal.inc'

  cosg=rnode(21,mptr)
  sing=rnode(22,mptr)

  k=MAX(kk,1)
  DO j=1,ny-1
    DO i=1,nx-1
      uc(i,j)=(cosg*(u(i,j,k)+u(i+1,j,k))-                              &
                 sing*(v(i,j,k)+v(i,j+1,k)))*0.5
      vc(i,j)=(sing*(u(i,j,k)+u(i+1,j,k))+                              &
                 cosg*(v(i,j,k)+v(i,j+1,k)))*0.5
    END DO
  END DO
  RETURN
END SUBROUTINE tranuv

  FUNCTION abssuv(array,i0,i1,j0,j1,k0,k1)
  REAL :: array(*)
  lmn=(i1-i0+1)*(j1-j0+1)*(k1-k0+1)
  abssuv=0.0
  DO ijk=1,lmn
    abssuv=abssuv+ABS(array(ijk))
  END DO
  RETURN
  END FUNCTION abssuv
