!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE UPDEXBC                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE updexbc(mptr,iexbcbuf)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate the external boundary data from base grid to grid mptr
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  04/15/1997
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!  mptr       Grid ID number
!
!  iexbcbuf   Pointer that points to the location of external
!             boundary data for base grid.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: lisvar(22)
  INTEGER :: mptr, iexbcbuf

  INTEGER :: nxs,nys,nzs, nxys, nxyzs
  INTEGER :: nx,ny,nz, nxy,nxyz
  INTEGER :: i,j,k,ijk

  INTEGER :: iu0exbc,iv0exbc,iw0exbc,ipt0exbc
  INTEGER :: iudtexbc,ivdtexbc,iwdtexbc,iptdtexbc
  INTEGER :: iuexbc,ivexbc,iwexbc,iptexbc

  INTEGER :: iptr1,iptr2,iptr3,iptr4,iptr5,iptr6

  INTEGER :: item3d,nwhts

  REAL :: delta_t
!
!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------
!
  INTEGER :: igetsp,igtexbc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( intrpodr == 2 ) THEN
    nwhts=9
  ELSE
    nwhts=4
  END IF

  nxs=node(5,1)
  nys=node(6,1)
  nzs=node(14,1)

  nx=node(5,mptr)
  ny=node(6,mptr)
  nz=node(14,mptr)

  nxy =nx*ny
  nxyz=nxy*nz

  nxys =nxs*nys
  nxyzs=nxys*nzs

  iptr1=igetsp(2*nxy)
  iptr2=igetsp(2*nxy)
  iptr3=igetsp(2*nxy)
  iptr4=igetsp(nxy)
  iptr5=igetsp(nwhts*nxy)
  iptr6=igetsp(2*nxy)

  item3d = igetsp( nxyzs )

  delta_t = curtim - ( abstfcst0 - abstinit )
!
!-----------------------------------------------------------------------
!
!  Interpolate uexbc
!
!-----------------------------------------------------------------------
!
  iu0exbc  = igtexbc(1,1,1)
  iudtexbc = igtexbc(1,12,1)

  DO k=1,nzs-1
    DO j=1,nys-1
      DO i=1,nxs
        ijk = (k-1)*nxs*nys + (j-1)*nxs + i
        a(item3d+ijk-1) = a(iu0exbc+ijk-1) + a(iudtexbc+ijk-1)*delta_t
      END DO
    END DO
  END DO

  iuexbc = iexbcbuf
  lisvar(1) = 1
  CALL itpnexbc(mptr,lisvar,1,item3d,iuexbc,                            &
       a(iptr1),a(iptr2),a(iptr3),                                      &
       a(iptr4),a(iptr5),a(iptr6),nxy,nwhts)
!
!-----------------------------------------------------------------------
!
!  Interpolate vexbc
!
!-----------------------------------------------------------------------
!
  iv0exbc  = igtexbc(1,2,1)
  ivdtexbc = igtexbc(1,13,1)

  DO k=1,nzs-1
    DO j=1,nys
      DO i=1,nxs-1
        ijk = (k-1)*nxs*nys + (j-1)*nxs + i
        a(item3d+ijk-1) = a(iv0exbc+ijk-1) + a(ivdtexbc+ijk-1)*delta_t
      END DO
    END DO
  END DO

  ivexbc = iexbcbuf + 1*nxyz
  lisvar(1) = 2
  CALL itpnexbc(mptr,lisvar,1,item3d,ivexbc,                            &
       a(iptr1),a(iptr2),a(iptr3),                                      &
       a(iptr4),a(iptr5),a(iptr6),nxy,nwhts)
!
!-----------------------------------------------------------------------
!
!  Interpolate wexbc
!
!-----------------------------------------------------------------------
!
  iw0exbc = igtexbc(1,3,1)
  iwdtexbc = igtexbc(1,14,1)

  DO k=1,nzs
    DO j=1,nys-1
      DO i=1,nxs-1
        ijk = (k-1)*nxs*nys + (j-1)*nxs + i
        a(item3d+ijk-1) = a(iw0exbc+ijk-1) + a(iwdtexbc+ijk-1)*delta_t
      END DO
    END DO
  END DO

  iwexbc = iexbcbuf + 2*nxyz
  lisvar(1) = 3
  CALL itpnexbc(mptr,lisvar,1,item3d,iwexbc,                            &
       a(iptr1),a(iptr2),a(iptr3),                                      &
       a(iptr4),a(iptr5),a(iptr6),nxy,nwhts)
!
!-----------------------------------------------------------------------
!
!  Interpolate scalar ptexbc.
!
!-----------------------------------------------------------------------
!
  ipt0exbc  = igtexbc(1,4,1)
  iptdtexbc = igtexbc(1,15,1)

  DO k=1,nzs-1
    DO j=1,nys-1
      DO i=1,nxs-1
        ijk = (k-1)*nxs*nys + (j-1)*nxs + i
        a(item3d+ijk-1) = a(ipt0exbc+ijk-1)+a(iptdtexbc+ijk-1)*delta_t
      END DO
    END DO
  END DO

  iptexbc = iexbcbuf + 3*nxyz
  lisvar(1) = 4

  CALL itpnexbc(mptr,lisvar,1,item3d,iptexbc,                           &
       a(iptr1),a(iptr2),a(iptr3),                                      &
       a(iptr4),a(iptr5),a(iptr6),nxy,nwhts)

  CALL retexbc(1,1,1,iu0exbc,   .true.)
  CALL retexbc(1,2,1,iv0exbc,   .true.)
  CALL retexbc(1,3,1,iw0exbc,   .true.)
  CALL retexbc(1,4,1,ipt0exbc,  .true.)
  CALL retexbc(1,12,1,iudtexbc, .true.)
  CALL retexbc(1,13,1,ivdtexbc, .true.)
  CALL retexbc(1,14,1,iwdtexbc, .true.)
  CALL retexbc(1,15,1,iptdtexbc,.true.)

  RETURN
END SUBROUTINE updexbc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE ITPNEXBC                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE itpnexbc(mrec,lisvar,nvar,isptr,irptr,                       &
           pts,irec,isrc,igsrc,whts,pts1,nptrec,nwhts)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate the external boundary data from base grid to grid mrec
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  04/15/1997
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: mrec,nptrec,nwhts
  INTEGER :: isptr,irptr, nvar
  INTEGER :: irec(2,nptrec),isrc(2,nptrec),igsrc(nptrec)
  INTEGER :: lisvar(nvar)

  REAL :: pts(2,nptrec),whts(nwhts,nptrec),pts1(2,nptrec)

  INTEGER :: nxr,nyr,nzr
  INTEGER :: nxs,nys,nzs
  INTEGER :: nx1,ny1
  INTEGER :: ist,jst,iend,jend

  INTEGER :: ivar
  INTEGER :: nw1,nw2
  INTEGER :: numpts,npts0
  INTEGER :: n3rd
  INTEGER :: ierfl,iptcmp

  INTEGER :: i,ii,j,jj, ip

  REAL :: cosr,sinr
  REAL :: x,y,xc,yc,xs,ys
  REAL :: xor,yor,xorp,yorp
  REAL :: dxr,dyr
  REAL :: xshift,yshift
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nw1=nint( SQRT(FLOAT(nwhts)) )
  nw2=nw1

  IF(nvar == 0) RETURN

  nxr=node(5,mrec)
  nyr=node(6,mrec)
  nzr=node(14,mrec)

  nxs=node(5,1)
  nys=node(6,1)
  nzs=node(14,1)

  cosr=rnode(21,mrec)
  sinr=rnode(22,mrec)
  dxr=rnode(9 ,mrec)
  dyr=rnode(10,mrec)
  xor=rnode(1,mrec)
  yor=rnode(2,mrec)

  ivar=lisvar(1)
  xshift=stgexbc(1,ivar)
  yshift=stgexbc(2,ivar)
  nx1=nxr-idmexbc(1,ivar)
  ny1=nyr-idmexbc(2,ivar)

  xc=xshift*dxr*cosr
  ys=yshift*dyr*sinr
  xs=xshift*dxr*sinr
  yc=yshift*dyr*cosr

  xorp=xor+xc-ys
  yorp=yor+xs+yc
!
! To calculate the absolute coordinates of all grid points including
! boundary and interior.
!
  ist=1
  jst=1
  iend=nx1
  jend=ny1
  DO j=jst,jend
    jj=(j-jst)*(iend-ist+1)
    DO i=ist,iend
      ii=jj+i-ist+1
      x=(i-1)*dxr
      y=(j-1)*dyr
      pts (1,ii) = xorp+x*cosr-y*sinr
      pts (2,ii) = yorp+x*sinr+y*cosr
      irec(1,ii) = i
      irec(2,ii) = j
    END DO
  END DO
  numpts=(iend-ist+1)*(jend-jst+1)
  npts0=numpts
!
! Find the source points for receive grid points pts
!
  ierfl=0
  iptcmp=1

  CALL getsrc(mrec,1, 1 ,xshift,yshift,pts,numpts,                      &
              irec,igsrc,ierfl,iptcmp)

  IF( (numpts /= npts0) .AND. verbose6 ) THEN
    WRITE(6,'('' WARNING: IN INTNSL,'',I4,'' POINTS ON GRID '',I4,      &
    &   '' got no source from source grid '',i4)')                      &
        npts0-numpts, mrec, 1
    WRITE(6,'('' IT WAS FOR EXBC VARIABLE NO. '',I4)') ivar
  END IF
!
! Calculate weights of interpolation for points pts.
!
  DO ip=1,numpts
    pts1(1,ip)=pts(1,ip)
    pts1(2,ip)=pts(2,ip)
  END DO

  CALL calwht_exbc(mrec,pts1,isrc,whts,numpts,nwhts,ivar)
!
! Interpolate for variables in list lisvar from source grid MSRC
! to receive grid MREC.
!
  DO i=1,nvar
    ivar=lisvar(i)
    n3rd=nzs

    CALL intrps(a(irptr),a(isptr),nxr,nyr,nxs,nys,n3rd,                 &
        irec,isrc,whts,numpts,nw1,nw2)

  END DO

  500   CONTINUE

  RETURN
END SUBROUTINE itpnexbc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE CALWHT_EXBC              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE calwht_exbc(mrec,points,isrc,whts,numpts,nwhts,ivar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  this routine calculates the weights for the interpolation
!  of mass points or vertical velocity (stuff that doesn't need
!  any rotation from base grid to grid mrec
!
!  if nw = 4 bilinear interp., = 9 bi-quadratic
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  04/15/1997
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!  mrec       Grid ID number
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: mrec
  INTEGER :: numpts,nwhts

  INTEGER :: isrc(2,numpts)
  REAL :: points(2,numpts),whts(nwhts,numpts)

  REAL :: xshift,yshift, yp
  REAL :: alpha,alpha2

  REAL :: tempx1,tempx2,tempx3
  REAL :: tempy1,tempy2,tempy3

  INTEGER :: i,ivar
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: cosc,sinc
  REAL :: dxc,dyc
  REAL :: odxc,odyc
  REAL :: xor,yor,xorp,yorp
  REAL :: xc,xs,yc,ys,x,y
!
!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------
!
  REAL :: fq, fl
  REAL :: xp,x1,x2,x3,x4

  fq(xp,x1,x2,x3,x4)=(xp-x2)*(xp-x3)/((x1-x2)*(x1-x3))+x4
  fl(xp,x1,x2)=(xp-x2)/(x1-x2)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  cosc = rnode( 21,1)
  sinc = rnode( 22,1)
  dxc  = rnode(9,1)
  dyc  = rnode(10,1)
  odxc = 1./dxc
  odyc = 1./dyc
  xor  = rnode(1,1)
  yor  = rnode(2,1)

  xshift=stgexbc(1,ivar)
  yshift=stgexbc(2,ivar)

  xc=xshift*dxc*cosc
  ys=yshift*dyc*sinc
  xs=xshift*dxc*sinc
  yc=yshift*dyc*cosc

  xorp=xor+xc-ys
  yorp=yor+xs+yc

  DO i=1,numpts
    x=points(1,i)
    y=points(2,i)
    points(1,i)=1+( cosc*(x-xorp)+sinc*(y-yorp))*odxc
    points(2,i)=1+(-sinc*(x-xorp)+cosc*(y-yorp))*odyc
  END DO

!
!   next, calculate weights for the interpolation
!
  IF ( nwhts == 9 ) THEN
    alpha=0.
    alpha2=0.
    DO i=1,numpts
      xp=points(1,i)-FLOAT( nint(points(1,i)) )
      yp=points(2,i)-FLOAT( nint(points(2,i)) )
      tempy1=fq(yp,-1.,0.,1.,alpha)
      tempy2=fq(yp,0.,-1.,1.,alpha2)
      tempy3=fq(yp,1.,-1.,0.,alpha)
      tempx1=fq(xp,-1.,0.,1.,alpha)
      tempx2=fq(xp,0.,-1.,1.,alpha2)
      tempx3=fq(xp,1.,-1.,0.,alpha)
      whts(1,i)=tempy1*tempx1
      whts(2,i)=tempy1*tempx2
      whts(3,i)=tempy1*tempx3
      whts(4,i)=tempy2*tempx1
      whts(5,i)=tempy2*tempx2
      whts(6,i)=tempy2*tempx3
      whts(7,i)=tempy3*tempx1
      whts(8,i)=tempy3*tempx2
      whts(9,i)=tempy3*tempx3
      isrc(1,i)=nint(points(1,i))
      isrc(2,i)=nint(points(2,i))
    END DO
  ELSE
    DO i=1,numpts
      xp=points(1,i)-FLOAT( INT(points(1,i)) )
      yp=points(2,i)-FLOAT( INT(points(2,i)) )
      tempy1=fl(yp,0.,1.)
      tempy2=fl(yp,1.,0.)
      tempx1=fl(xp,0.,1.)
      tempx2=fl(xp,1.,0.)
      whts(1,i)=tempy1*tempx1
      whts(2,i)=tempy1*tempx2
      whts(3,i)=tempy2*tempx1
      whts(4,i)=tempy2*tempx2
      isrc(1,i)=INT(points(1,i))
      isrc(2,i)=INT(points(2,i))
    END DO
  END IF

  RETURN
END SUBROUTINE calwht_exbc
