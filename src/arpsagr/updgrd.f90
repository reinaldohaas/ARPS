!
!--------------------------------------------------------------------------
!

SUBROUTINE updgrd(mptr)
!
! To update a coarse grid by interpolating from finer grid(s).
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
  DIMENSION msrcs(30),iupvar(nxy2d+nxyz3d),iupsta(nxy2d+nxyz3d)
  LOGICAL :: smstg2

  IF( intrpodr == 2 ) THEN
    nwhts=9
  ELSE
    nwhts=4
  END IF

  level = node(4,mptr)
  IF(level >= lfine) RETURN
!
! Get all the possible source grids on the next level for mptr.
!
  numscs=0
  mgrid=lstart(level+1)
  15    IF(mgrid == 0) GO TO 20
  numscs=numscs+1
  msrcs(numscs)=mgrid
  mgrid=node(10,mgrid)
  GO TO 15
  20    CONTINUE

  nx=node(5,mptr)
  ny=node(6,mptr)
!
! --------------------------------------------------------------------
! update 2-d scalars
! --------------------------------------------------------------------
!
  IF (verbose6) WRITE(6,'(''  UPDATING 2-D SCALAR(S):''//)')
  ndim=2
!
! initialize the updating status flag with value 1
!
  DO i=1,nxy2d
    iupsta(i)=1
  END DO

  100   CONTINUE
!
! search for unupdated 2-d scalars and group them according to staggering
!
  nupvar=0
  DO iup=1,nxy2d
!
! look for 1st 2-d scalar that needs updating
!
    IF(iupxy(1,iup) == 1.AND.iupsta(iup) == 1.AND. iupxy(2,iup) == 0)THEN
      nupvar=nupvar+1
      iupvar(nupvar)=iup
      iupsta(iup)=0
      DO ivar=iup+1,nxy2d
!
! list additonal scalars with matching staggering to that of iup.
!
        IF(iupxy(1,ivar) == 1.AND.iupsta(ivar) == 1.AND.                &
              iupxy(2,ivar) == 0.AND. smstg2(ivar,iup,2)) THEN
          nupvar=nupvar+1
          iupvar(nupvar)=ivar
          iupsta(ivar)=0
        END IF
      END DO
    END IF
  END DO
!
! if no 2-d scalar needs updating goto next step
!
  IF(nupvar == 0)GO TO 200
!
! perform interpolations for 3-d scalars listed in iupvar(nupvar).
!
  IF( verbose6) THEN
    WRITE(6,'(/''  UPDATING '',I3,''  2-D SCALAR(S)''/)')nupvar
    WRITE(6,9992) (iupvar(ii),ii=1,nupvar)
9992 FORMAT(1X, 10I5)
  END IF
!
  mxpts=nx*ny
  iptr1=igetsp(2*mxpts)
  iptr2=igetsp(2*mxpts)
  iptr3=igetsp(2*mxpts)
  iptr4=igetsp(  mxpts)
  iptr5=igetsp(nwhts*mxpts)
  iptr6=igetsp(2*mxpts)
  iptr7=igetsp(2*mxpts)
!
  CALL itpscl(mptr,msrcs,numscs,a(iptr1),a(iptr2),a(iptr3),             &
              a(iptr4),a(iptr5),a(iptr6),a(iptr7),                      &
              mxpts,nwhts,iupvar,nupvar,ndim)
!
!  restore tem space
!
  CALL resett
!
! go back to check if any more var's need updating
!
  GO TO 100
  200   CONTINUE
!
! --------------------------------------------------------------------
! update 3-d scalars
! --------------------------------------------------------------------
!
  IF(verbose6) WRITE(6,'(''  UPDATING 3-D SCALAR(S):''//)')
  ndim=3
!
! initialize the updating status flag with value 1
!
  DO i=1,nxyz3d
    iupsta(i)=1
  END DO
!
! search for unupdated 3-d scalars and group them according to staggering
!
  300   CONTINUE
  nupvar=0
  DO iup=1,nxyz3d
!
! look for 1st 3-d scalar that needs updating
!
    IF(iupxyz(1,iup) == 1.AND.iupsta(iup) == 1.AND. iupxyz(2,iup) == 0)THEN
      nupvar=nupvar+1
      iupvar(nupvar)=iup
      iupsta(iup)=0
      DO ivar=iup+1,nxyz3d
!
! list additonal 3-d scalars with matching staggering to that of iup.
!
        IF(iupxyz(1,ivar) == 1.AND.iupsta(ivar) == 1.AND.               &
              iupxyz(2,ivar) == 0.AND. smstg2(ivar,iup,3))THEN
          nupvar=nupvar+1
          iupvar(nupvar)=ivar
          iupsta(ivar)=0
        END IF
      END DO
    END IF
  END DO
!
! if no more 3-d scalar needs updating goto next step
!
  IF(nupvar == 0) GO TO 400
!
! perform interpolations for 3-d scalars listed in iupvar(nupvar).
!
  IF(verbose6) THEN
    WRITE(6,'(/''  UPDATING '',I3,''  3-D SCALAR(S)''/)')nupvar
    WRITE(6,9992) (iupvar(ii),ii=1,nupvar)
  END IF
!
  mxpts=nx*ny
  iptr1=igetsp(2*mxpts)
  iptr2=igetsp(2*mxpts)
  iptr3=igetsp(2*mxpts)
  iptr4=igetsp(  mxpts)
  iptr5=igetsp(nwhts*mxpts)
  iptr6=igetsp(2*mxpts)
  iptr7=igetsp(2*mxpts)
!
  CALL itpscl(mptr,msrcs,numscs,a(iptr1),a(iptr2),a(iptr3),             &
              a(iptr4),a(iptr5),a(iptr6),a(iptr7),                      &
              mxpts,nwhts,iupvar,nupvar,ndim)
!
!  restore tem space
!
  CALL resett
!
! go back to check if any more var's need updating
!
  GO TO 300
  400   CONTINUE
!
!-------------------------------------------------------------------
! update 2-d vectors
! --------------------------------------------------------------------
!
!  write(6,'(''  updating 2-d vector(s):''//)')
  ndim=2
!
! initialize the updating status flag with value 1
!
  DO i=1,nxy2d
    iupsta(i)=1
  END DO

  500   CONTINUE
!
! search for a vector that needs updating and store the compnt No. in
! array iupvar.
!
  nupvar=0
  DO iup=1,nxy2d
!
! look for 1st 2-d vector that needs updating
!
    IF(iupxy(1,iup) == 1.AND.iupxy(2,iup) > 0.AND. iupsta(iup) == 1) THEN
      iup2=iupxy(2,iup)
      iupvar(nupvar+1)=iup
      iupvar(nupvar+2)=iup2
      nupvar=nupvar+2
      iupsta(iup )=0
      iupsta(iup2)=0
      EXIT
    END IF
  END DO
  67    CONTINUE
!
! if no more 2-d vector needs updating goto next step
!
  IF(nupvar == 0)GO TO 600
!
! perform interpolations for 2-d vectors listed in iupvar(nupvar).
!
  IF(verbose6) THEN
    WRITE(6,'(//''  UPDATING 2-D VECTOR '',I5,I5//)')                   &
         iupvar(1),iupvar(2)
  END IF
!
  mxpts=nx*ny
  iptr1=igetsp(4*mxpts)
  iptr2=igetsp(4*mxpts)
  iptr3=igetsp(4*mxpts)
  iptr4=igetsp(2*mxpts)
  iptr5=igetsp(2*nwhts*mxpts)
  iptr6=igetsp(4*mxpts)
  iptr7=igetsp(4*mxpts)
!
  CALL itpvtr(mptr,msrcs,numscs,a(iptr1),a(iptr2),a(iptr3),             &
              a(iptr4),a(iptr5),a(iptr6),a(iptr7),                      &
              mxpts,nwhts,iupvar,nupvar,ndim)
!
!  restore tem space
!
  CALL resett
!
! go back to check if any more var's need updating
!
  GO TO 500
  600   CONTINUE
!
!------------------------------------------------------------------
! update 3-d vectors
!------------------------------------------------------------------
!
  IF(verbose6) WRITE(6,'(''  UPDATING 3-D VECTOR(S):''//)')
  ndim=3
!
! initialize the updating status flag with value 1
!
  DO i=1,nxyz3d
    iupsta(i)=1
  END DO
!
  700   CONTINUE
!
! search for a vector that needs updating and store the compnt No. in
! array iupvar.
!
  nupvar=0
  DO iup=1,nxyz3d
!
! pick out the 1st 3-d vector that needs updating
!
    IF(iupxyz(1,iup) == 1.AND.iupxyz(2,iup) > 0.AND. iupsta(iup) == 1) THEN
      iup2=iupxyz(2,iup)
      iupvar(nupvar+1)=iup
      iupvar(nupvar+2)=iup2
      nupvar=nupvar+2
      iupsta(iup )=0
      iupsta(iup2)=0
      EXIT
    END IF
  END DO
  92    CONTINUE
!
! if no more 3-d vector needs updating goto next step
!
  IF(nupvar == 0)GO TO 800
!
! perform interpolations for 3-d vectors listed in iupvar(nupvar).
!
  IF(verbose6) THEN
    WRITE(6,'(//''  UPDATING 3-D VECTOR '',I5,I5//)')                   &
          iupvar(1),iupvar(2)
  END IF
!
  mxpts=nx*ny
  iptr1=igetsp(4*mxpts)
  iptr2=igetsp(4*mxpts)
  iptr3=igetsp(4*mxpts)
  iptr4=igetsp(2*mxpts)
  iptr5=igetsp(2*nwhts*mxpts)
  iptr6=igetsp(4*mxpts)
  iptr7=igetsp(4*mxpts)
!
  CALL itpvtr(mptr,msrcs,numscs,a(iptr1),a(iptr2),a(iptr3),             &
              a(iptr4),a(iptr5),a(iptr6),a(iptr7),                      &
              mxpts,nwhts,iupvar,nupvar,ndim)
!
!  restore tem space
!
  CALL resett
!
! go back to check if any more var's need updating
!
  GO TO 700
  800   CONTINUE
  RETURN
END SUBROUTINE updgrd
!
!--------------------------------------------------------------------------
!

SUBROUTINE itpscl(mrec,msrcs,numscs,ptsrec,                             &
           irec,isrc,igsrc,whts,pts,irecp,mxpts,nwhts,                  &
           iupvar,nupvar,ndim)
!
  DIMENSION msrcs(numscs),ptsrec(2,mxpts),irec(2,mxpts),                &
            isrc(2,mxpts),igsrc(mxpts),whts(nwhts,mxpts),               &
            pts(2,mxpts),irecp(2,mxpts),                                &
            iupvar(nupvar)

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  nw1=nint( SQRT(FLOAT(nwhts)) )
  nw2=nw1
!
! Calculate and store the coor's of p-points of the receive
! grid in array ptsrec with their indecies stored in irec.
!
  IF(nupvar == 0) RETURN

  nx=node(5,mrec)
  ny=node(6,mrec)
  cosr=rnode(21,mrec)
  sinr=rnode(22,mrec)
  dxr=rnode(9 ,mrec)
  dyr=rnode(10,mrec)
  xor=rnode(1,mrec)
  yor=rnode(2,mrec)

  ivar=iupvar(1)
  IF(ndim == 2) THEN
    xshift=stgxy(1,ivar)
    yshift=stgxy(2,ivar)
    nx1=nx-idmxy(1,ivar)
    ny1=ny-idmxy(2,ivar)
  ELSE IF(ndim == 3)THEN
    xshift=stgxyz(1,ivar)
    yshift=stgxyz(2,ivar)
    nx1=nx-idmxyz(1,ivar)
    ny1=ny-idmxyz(2,ivar)
  ELSE
    PRINT*,'Error in ITPSCL: Variable dimension ndim set wrong!'
    PRINT*,'ndim should be 2 or 3, input was',ndim
  END IF

  xc=xshift*dxr*cosr
  ys=yshift*dyr*sinr
  xs=xshift*dxr*sinr
  yc=yshift*dyr*cosr

  xorp=xor+xc-ys
  yorp=yor+xs+yc
!
! To calculate the absolute coordinates of the model interior
! points that may need updating.
!
! Attention:
!
! here we assume that in x direction the interior points start
! at i=2 and finish at i=nx1-1. The boundaries are at i=1 and nx1.
! Similarly in y direction, the boundaries are at j=1 and ny1.
! If in a particular solver (model), the boundaries for certain
! variables are defined at points other the above, the use needs to
! change the following loop accordingly (i.e. change ist, jst, iend
! and jend). As a result, the part in routine ITPBSL that calculates
! the points for boundary value interpolations also need change.
!
  ist=2
  jst=2
  iend=nx1-1
  jend=ny1-1
  DO j=jst,jend
    jj=(j-jst)*(iend-ist+1)
    DO i=ist,iend
      ii           = jj+i-ist+1
      x            = (i-1)*dxr
      y            = (j-1)*dyr
      ptsrec(1,ii) = xorp+x*cosr-y*sinr
      ptsrec(2,ii) = yorp+x*sinr+y*cosr
      irec  (1,ii) = i
      irec  (2,ii) = j
    END DO
  END DO
  numpts=(iend-ist+1)*(jend-jst+1)
!
! Find the source points for receive grid points ptsrec
!
  ierfl=0
  iptcmp=1
  CALL getsrc(mrec,msrcs,numscs,xshift,yshift,ptsrec,numpts,            &
       irec,igsrc,ierfl,iptcmp)
!
! if no point can be updated from listed source grids, skip through
! this loop
!
  IF(numpts == 0) THEN
    WRITE(6,'(''NO SRC. POINT FOUND FOR'',I2,''-d variable(s)'')')ndim
    WRITE(6,'(20i4)') iupvar
    RETURN
  END IF
!
! Loop through all source grids
!
  DO ks=1,numscs

    msrc=msrcs(ks)
!
! sort rec. grid points according to their source grid
!
    npts=0
    DO np=1,numpts
      IF(igsrc(np) == msrc) THEN
        npts=npts+1
        pts(1,npts)=ptsrec(1,np)
        pts(2,npts)=ptsrec(2,np)
        irecp(1,npts)=irec(1,np)
        irecp(2,npts)=irec(2,np)
      END IF
    END DO
    IF(npts == 0) CYCLE

    nxr=node(5,mrec)
    nyr=node(6,mrec)
    nzr=node(14,mrec)

    nxs=node(5,msrc)
    nys=node(6,msrc)
    nzs=node(14,msrc)

!
! Calculate weights of interpolation for points pts.
!
    CALL calwht(mrec,msrc,pts,isrc,whts,npts,nwhts,ivar,ndim)
!  print*,'after calling calwht'
!    do 1111 ii=1,npts
!     write(6,'(5i5,10f10.5)') ii,irecp(1,ii),irecp(2,ii),
!    : isrc(1,ii),isrc(2,ii),(whts(jj,ii),jj=1,nw1*nw2)
!1111   continue
!
! update variables in list iupvar
!
    DO i=1,nupvar
      ivar=iupvar(i)
      IF(ndim == 2) THEN
        n3rd=1
        irptr=igtnxy(mrec,ivar,1)
        isptr=igtnxy(msrc,ivar,1)
      ELSE
        n3rd=nzs
        irptr=igtxyz(mrec,ivar,1)
        isptr=igtxyz(msrc,ivar,1)
      END IF
      needsp =nxs*nys*n3rd
      iavrg=igetsp(needsp)

      CALL hrzavg(a(iavrg),a(isptr),nxs,nys,n3rd)

      CALL intrps(a(irptr),a(iavrg),nxr,nyr,nxs,nys,n3rd,               &
          irecp,isrc,whts,npts,nw1,nw2)

      IF(ndim == 2) THEN
        CALL retnxy(mrec,ivar,1,irptr,.true.)
        CALL retnxy(msrc,ivar,1,isptr,.false.)
      ELSE
        CALL retxyz(mrec,ivar,1,irptr,.true.)
        CALL retxyz(msrc,ivar,1,isptr,.false.)
      END IF
      CALL reclam(iavrg,nxs*nys*n3rd)
    END DO
!
  END DO
  500   CONTINUE
  RETURN
END SUBROUTINE itpscl
!
!--------------------------------------------------------------------------
!

SUBROUTINE intrps(varrec,varsrc,nxr,nyr,nxs,nys,nz,                     &
           irec,isrc,whts,numpts,nw1,nw2)
  DIMENSION varrec(nxr,nyr,nz),varsrc(nxs,nys,nz),                      &
            irec(2,numpts),isrc(2,numpts),whts(nw1,nw2,numpts)
!
  itest = 0
  IF(itest == 1) THEN
    WRITE(6,'(''  DATA DUMP IN INTRPS '',8I6)')                         &
        nxr,nyr,nxs,nys,nz,numpts,nw1,nw2
    DO ip=1,20
      WRITE(7,1114) irec(1,ip),irec(2,ip),isrc(1,ip),isrc(2,ip),        &
            whts(1,1,ip),whts(2,1,ip),whts(1,2,ip),whts(2,2,ip)
    END DO
1114 FORMAT(1X,4I7,4E10.3)
  END IF
  IF(nw1 == 2.AND.nw2 == 2) THEN
!
! The interpolation calculations are all independent of each other,
! they can be vecotrized in inner loop and parellelized in outer loop.
!
    nzswap=2
!
! if nz if small, swap the loop order
!
    IF(nz <= nzswap)THEN
!fpp$ nodepchk
      DO ks=1,nz
!DIR$ IVDEP
        DO np=1,numpts
          ir = irec(1,np)
          jr = irec(2,np)
          is = isrc(1,np)
          js = isrc(2,np)
          varrec(ir,jr,ks)=varsrc(is,js,ks)*whts(1,1,np)                &
              +varsrc(is+1,js,ks)*whts(2,1,np)                          &
              +varsrc(is,js+1,ks)*whts(1,2,np)                          &
              +varsrc(is+1,js+1,ks)*whts(2,2,np)
        END DO
      END DO
    ELSE
!fpp$ nodepchk
      DO np=1,numpts
        ir = irec(1,np)
        jr = irec(2,np)
        is = isrc(1,np)
        js = isrc(2,np)
!DIR$ IVDEP
        DO ks=1,nz
          varrec(ir,jr,ks)=varsrc(is,js,ks)*whts(1,1,np)                &
              +varsrc(is+1,js,ks)*whts(2,1,np)                          &
              +varsrc(is,js+1,ks)*whts(1,2,np)                          &
              +varsrc(is+1,js+1,ks)*whts(2,2,np)
        END DO
      END DO
    END IF

  ELSE

    IF(nz <= nzswap)THEN
!fpp$ nodepchk
      DO ks=1,nz
!DIR$ IVDEP
        DO np=1,numpts
          ir = irec(1,np)
          jr = irec(2,np)
          is = isrc(1,np)
          js = isrc(2,np)
          varrec(ir,jr,ks)=                                             &
              varsrc(is-1,js-1,ks)*whts(1,1,np)                         &
              +varsrc(is  ,js-1,ks)*whts(2,1,np)                        &
              +varsrc(is+1,js-1,ks)*whts(3,1,np)                        &
              +varsrc(is-1,js  ,ks)*whts(1,2,np)                        &
              +varsrc(is  ,js  ,ks)*whts(2,2,np)                        &
              +varsrc(is+1,js  ,ks)*whts(3,2,np)                        &
              +varsrc(is-1,js+1,ks)*whts(1,3,np)                        &
              +varsrc(is  ,js+1,ks)*whts(2,3,np)                        &
              +varsrc(is+1,js+1,ks)*whts(3,3,np)
        END DO
      END DO
    ELSE
!fpp$ nodepchk
      DO np=1,numpts
        ir = irec(1,np)
        jr = irec(2,np)
        is = isrc(1,np)
        js = isrc(2,np)
!DIR$ IVDEP
        DO ks=1,nz
          varrec(ir,jr,ks)=                                             &
              varsrc(is-1,js-1,ks)*whts(1,1,np)                         &
              +varsrc(is  ,js-1,ks)*whts(2,1,np)                        &
              +varsrc(is+1,js-1,ks)*whts(3,1,np)                        &
              +varsrc(is-1,js  ,ks)*whts(1,2,np)                        &
              +varsrc(is  ,js  ,ks)*whts(2,2,np)                        &
              +varsrc(is+1,js  ,ks)*whts(3,2,np)                        &
              +varsrc(is-1,js+1,ks)*whts(1,3,np)                        &
              +varsrc(is  ,js+1,ks)*whts(2,3,np)                        &
              +varsrc(is+1,js+1,ks)*whts(3,3,np)
        END DO
      END DO
    END IF

  END IF

  RETURN
END SUBROUTINE intrps
!
!--------------------------------------------------------------------------
!

SUBROUTINE itpvtr(mrec,msrcs,numscs,ptsrec,                             &
           irec,isrc,igsrc,whts,pts,irecp,mxpts,nwhts,                  &
           iupvar,nupvar,ndim)
!
  DIMENSION msrcs(numscs),ptsrec(2,mxpts,2),irec(2,mxpts,2),            &
            isrc(2,mxpts,2),igsrc(mxpts,2),whts(nwhts,mxpts,2),         &
            pts(2,mxpts,2),irecp(2,mxpts,2),iupvar(nupvar),             &
            numpts(2),npts(2)

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  LOGICAL :: samstg,smstg2

  nw1=nint( SQRT(FLOAT(nwhts)) )
  nw2=nw1
!
  IF(nupvar == 0) RETURN
  nxr=node(5,mrec)
  nyr=node(6,mrec)
  nzr=node(14,mrec)
  cosr=rnode(21,mrec)
  sinr=rnode(22,mrec)
  dxr=rnode(9 ,mrec)
  dyr=rnode(10,mrec)
  xor=rnode(1,mrec)
  yor=rnode(2,mrec)
!
  CALL chkdim(ndim,'ITPVTR')
!
! To loop through all vectors, which has two components stored consecutively
! in nupvar. The 1st compt. has to be the one in x direction.
!
  DO ivtr=1,nupvar,2

    ivar1=iupvar(ivtr)
    ivar2=iupvar(ivtr+1)

    samstg=smstg2(ivar1,ivar2,ndim)
    IF(ndim == 2) THEN
      xshft = MAX( stgxy(1,ivar1),stgxy(1,ivar2) )
      yshft = MAX( stgxy(2,ivar1),stgxy(2,ivar2) )
    ELSE IF(ndim == 3)THEN
      xshft = MAX( stgxyz(1,ivar1),stgxyz(1,ivar2) )
      yshft = MAX( stgxyz(2,ivar1),stgxyz(2,ivar2) )
    END IF
!
! Find the points on rec. grid for both compnts. and the source grid for
! them.If the two compnts have the samew staggering, they have the
! same coordinate points and the source grid, therefore those for the
! 2nd compnt. do not need to be recalculated.
!
    DO icmpr=1,2
!
! icmpr=1 for the 1st component of vectors
!   =2 for the 2nd component of vectors
!
      IF(icmpr == 2.AND.samstg)THEN
        numpts(2)=numpts(1)
        DO ip=1,numpts(2)
          ptsrec(1,ip,2)=ptsrec(1,ip,1)
          ptsrec(2,ip,2)=ptsrec(2,ip,1)
          irec(1,ip,2)=irec(1,ip,1)
          irec(2,ip,2)=irec(2,ip,1)
          igsrc(ip,2)=igsrc(ip,1)
        END DO
        CYCLE
      END IF

      ivar=iupvar(ivtr+icmpr-1)

      IF(ndim == 2) THEN
        xshift=stgxy(1,ivar)
        yshift=stgxy(2,ivar)
        nx1=nxr-idmxy(1,ivar)
        ny1=nyr-idmxy(2,ivar)
      ELSE IF(ndim == 3)THEN
        xshift=stgxyz(1,ivar)
        yshift=stgxyz(2,ivar)
        nx1=nxr-idmxyz(1,ivar)
        ny1=nyr-idmxyz(2,ivar)
      END IF

      xc=xshift*dxr*cosr
      ys=yshift*dyr*sinr
      xs=xshift*dxr*sinr
      yc=yshift*dyr*cosr

      xorp=xor+xc-ys
      yorp=yor+xs+yc
!
! The following loop calculates the absolute coordinates of the model
! interior points that may need updating.
!
! Attention:
!
! here we assume that in x direction the interior points start
! at i=2 and finish at i=nx1-1. The boundaries are at i=1 and nx1.
! Similarly in y direction, the boundaries are at j=1 and ny1.
! If in a particular solver (model), the boundaries for certain
! variables are defined at points other the above, the use needs to
! change the following loop accordingly (i.e. change ist, jst, iend
! and jend). As a result, the part in routine ITPBVT that calculates
! the points for boundary value interpolations also need change.
!
      ist=2
      jst=2
      iend=nx1-1
      jend=ny1-1
      DO j=jst,jend
        jj=(j-jst)*(iend-ist+1)
        DO i=ist,iend
          ii=jj+i-ist+1
          x=(i-1)*dxr
          y=(j-1)*dyr
          ptsrec(1,ii,icmpr) = xorp+x*cosr-y*sinr
          ptsrec(2,ii,icmpr) = yorp+x*sinr+y*cosr
          irec  (1,ii,icmpr) = i
          irec  (2,ii,icmpr) = j
        END DO
      END DO
      numpts(icmpr)=(iend-ist+1)*(jend-jst+1)
      16    CONTINUE
!
! Find the source grid for receive grid points ptsrec
!
      ierfl=0
      iptcmp=1
      CALL getsrc(mrec,msrcs,numscs,xshft,yshft,ptsrec(1,1,icmpr),      &
          numpts(icmpr),irec(1,1,icmpr),igsrc(1,icmpr),ierfl,iptcmp)

    END DO
!
! if no point can be updated from listed source grids, skip through
! this loop
!
    IF(numpts(1) == 0.AND.numpts(2) == 0) THEN
      WRITE(6,'(''NO SRC. POINT FOUND FOR'',I2,                         &
      &   ''-d vector no. '', 2I6)')ndim,iupvar(ivtr),iupvar(ivtr+1)
      CYCLE
    END IF
!
! find storage pointers for the receive grid variables (2 compnts of vector).
! uppack them if neccesary.
!
    IF(ndim == 2) THEN
      irptr1=igtnxy(mrec,ivar1,1)
      irptr2=igtnxy(mrec,ivar2,1)
    ELSE
      irptr1=igtxyz(mrec,ivar1,1)
      irptr2=igtxyz(mrec,ivar2,1)
    END IF
!
! Loop through all source grids
!
    DO ks=1,numscs

      msrc=msrcs(ks)
      nxs=node(5,msrc)
      nys=node(6,msrc)
      nzs=node(14,msrc)
      coss=rnode(21,msrc)
      sins=rnode(22,msrc)

      cossr= coss*cosr+sins*sinr
      sinsr= coss*sinr-sins*cosr

      DO icmps=1,2

        DO  icmpr=1,2

          IF(icmpr == 2.AND.samstg)THEN
            npts(2)=npts(1)
            DO ip=1,npts(2)
              pts  (1,ip,2)=pts  (1,ip,1)
              pts  (2,ip,2)=pts  (2,ip,1)
              isrc (1,ip,2)=isrc (1,ip,1)
              isrc (2,ip,2)=isrc (2,ip,1)
              irecp(1,ip,2)=irecp(1,ip,1)
              irecp(2,ip,2)=irecp(2,ip,1)
            END DO
            DO iw=1,nwhts
              DO ip=1,npts(2)
                whts(iw,ip,2)=whts(iw,ip,1)
              END DO
            END DO
            PRINT*,'irecp, isrc, whts with icmps=:',icmps
            CYCLE
          END IF
!
! sort rec. grid points according to their source grid
!
          npt=0
          DO np=1,numpts(icmpr)
            IF(igsrc(np,icmpr) == msrc) THEN
              npt=npt+1
              pts(1,npt,icmpr)=ptsrec(1,np,icmpr)
              pts(2,npt,icmpr)=ptsrec(2,np,icmpr)
              irecp(1,npt,icmpr)=irec(1,np,icmpr)
              irecp(2,npt,icmpr)=irec(2,np,icmpr)
            END IF
          END DO
          npts(icmpr)=npt
          IF(npts(icmpr) == 0) CYCLE
!
! Calculate weights of interpolation for points pts.
!
          ivsrc=iupvar(ivtr+icmps-1)
          CALL calwht(mrec,msrc,pts(1,1,icmpr),isrc(1,1,icmpr),         &
                      whts(1,1,icmpr),npts(icmpr),nwhts,ivsrc,ndim)
        END DO
        IF(npts(1) == 0.AND.npts(2) == 0) CYCLE
!
! find pointer for one of the components (icmps) of the source grid
! vector. If it is packed, do uppacking for this variable.
!
        ivsrc=iupvar(ivtr+icmps-1)
        IF(ndim == 2) THEN
          n3rd=1
          isptr=igtnxy(msrc,ivsrc,1)
        ELSE
          n3rd=nzs
          isptr=igtxyz(msrc,ivsrc,1)
        END IF
        needsp = nxs*nys*n3rd
        iavrg=igetsp(needsp)
        CALL hrzavg(a(iavrg),a(isptr),nxs,nys,n3rd)

        IF(icmps == 1)THEN
          proj1=cossr
          proj2=-sinsr
        ELSE
          proj1=sinsr
          proj2=cossr
        END IF

        icmpr=1
        IF(npts(icmpr) /= 0) THEN
          CALL intrpv(a(irptr1),a(iavrg),nxr,nyr,nxs,nys,n3rd,          &
              irecp(1,1,icmpr),isrc(1,1,icmpr),whts(1,1,icmpr),         &
              npts(icmpr),nw1,nw2,icmps,proj1)
        END IF

        icmpr=2
        IF(npts(icmpr) /= 0) THEN
          CALL intrpv(a(irptr2),a(iavrg),nxr,nyr,nxs,nys,n3rd,          &
              irecp(1,1,icmpr),isrc(1,1,icmpr),whts(1,1,icmpr),         &
              npts(icmpr),nw1,nw2,icmps,proj2)
        END IF

        CALL reclam(iavrg,nxs*nys*n3rd)
        IF(ndim == 2) THEN
          CALL retnxy(msrc,ivsrc,1,isptr,.false.)
        ELSE
          CALL retxyz(msrc,ivsrc,1,isptr,.false.)
        END IF
      END DO
!
    END DO
    IF(ndim == 2) THEN
      CALL retnxy(mrec,ivar1,1,irptr1,.true.)
      CALL retnxy(mrec,ivar2,1,irptr2,.true.)
    ELSE
      CALL retxyz(mrec,ivar1,1,irptr1,.true.)
      CALL retxyz(mrec,ivar2,1,irptr2,.true.)
    END IF
  END DO
  RETURN
END SUBROUTINE itpvtr
!
!--------------------------------------------------------------------------
!

SUBROUTINE intrpv(varrec,varsrc,nxr,nyr,nxs,nys,nz,                     &
           irec,isrc,whts,numpts,nw1,nw2,mode,proj)
  DIMENSION varrec(nxr,nyr,nz),varsrc(nxs,nys,nz),                      &
            irec(2,numpts),isrc(2,numpts),whts(nw1,nw2,numpts)
!
  itest = 0
  IF(itest == 1) THEN
    WRITE(6,'(''  DATA DUMP IN INTRPV '',8I6)')                         &
        nxr,nyr,nxs,nys,nz,numpts,nw1,nw2
    DO ip=1,20
      WRITE(7,1114) irec(1,ip),irec(2,ip),isrc(1,ip),isrc(2,ip),        &
            whts(1,1,ip),whts(2,1,ip),whts(1,2,ip),whts(2,2,ip)
    END DO
1114 FORMAT(1X,4I7,4E10.3)
  END IF

  iaccum=0
  IF(mode == 2) iaccum=1

  IF(nw1 == 2.AND.nw2 == 2) THEN
!
! The interpolation calculations are all independent of each other,
! they can be vecotrized in inner loop and parellelized in outer loop.
!
    nzswap=2
!
! if nz if small, swap the loop order
!
    IF(nz <= nzswap)THEN
!fpp$ nodepchk
      DO ks=1,nz
!DIR$ IVDEP
        DO np=1,numpts
          ir = irec(1,np)
          jr = irec(2,np)
          is = isrc(1,np)
          js = isrc(2,np)
          varrec(ir,jr,ks)=varrec(ir,jr,ks)*iaccum                      &
              +(varsrc(is,js,ks)*whts(1,1,np)                           &
              +varsrc(is+1,js,ks)*whts(2,1,np)                          &
              +varsrc(is,js+1,ks)*whts(1,2,np)                          &
              +varsrc(is+1,js+1,ks)*whts(2,2,np))*proj
        END DO
      END DO
    ELSE
!fpp$ nodepchk
      DO np=1,numpts
        ir = irec(1,np)
        jr = irec(2,np)
        is = isrc(1,np)
        js = isrc(2,np)
!DIR$ IVDEP
        DO ks=1,nz
          varrec(ir,jr,ks)=varrec(ir,jr,ks)*iaccum                      &
              +(varsrc(is,js,ks)*whts(1,1,np)                           &
              +varsrc(is+1,js,ks)*whts(2,1,np)                          &
              +varsrc(is,js+1,ks)*whts(1,2,np)                          &
              +varsrc(is+1,js+1,ks)*whts(2,2,np))*proj
        END DO
      END DO
    END IF

  ELSE

    IF(nz <= nzswap)THEN
!fpp$ nodepchk
      DO ks=1,nz
!DIR$ IVDEP
        DO np=1,numpts
          ir = irec(1,np)
          jr = irec(2,np)
          is = isrc(1,np)
          js = isrc(2,np)
          varrec(ir,jr,ks)=varrec(ir,jr,ks)*iaccum                      &
              +(varsrc(is-1,js-1,ks)*whts(1,1,np)                       &
              +varsrc(is  ,js-1,ks)*whts(2,1,np)                        &
              +varsrc(is+1,js-1,ks)*whts(3,1,np)                        &
              +varsrc(is-1,js  ,ks)*whts(1,2,np)                        &
              +varsrc(is  ,js  ,ks)*whts(2,2,np)                        &
              +varsrc(is+1,js  ,ks)*whts(3,2,np)                        &
              +varsrc(is-1,js+1,ks)*whts(1,3,np)                        &
              +varsrc(is  ,js+1,ks)*whts(2,3,np)                        &
              +varsrc(is+1,js+1,ks)*whts(3,3,np))*proj
        END DO
      END DO
    ELSE
!fpp$ nodepchk
      DO np=1,numpts
        ir = irec(1,np)
        jr = irec(2,np)
        is = isrc(1,np)
        js = isrc(2,np)
!DIR$ IVDEP
        DO ks=1,nz
          varrec(ir,jr,ks)=varrec(ir,jr,ks)*iaccum                      &
              +(varsrc(is-1,js-1,ks)*whts(1,1,np)                       &
              +varsrc(is  ,js-1,ks)*whts(2,1,np)                        &
              +varsrc(is+1,js-1,ks)*whts(3,1,np)                        &
              +varsrc(is-1,js  ,ks)*whts(1,2,np)                        &
              +varsrc(is  ,js  ,ks)*whts(2,2,np)                        &
              +varsrc(is+1,js  ,ks)*whts(3,2,np)                        &
              +varsrc(is-1,js+1,ks)*whts(1,3,np)                        &
              +varsrc(is  ,js+1,ks)*whts(2,3,np)                        &
              +varsrc(is+1,js+1,ks)*whts(3,3,np))*proj
        END DO
      END DO
    END IF

  END IF

  RETURN
END SUBROUTINE intrpv
!
!--------------------------------------------------------------------------
!

SUBROUTINE getsrc(mrec,msrcs,numscs,xshift,yshift,ptsrec,numpts,        &
           irec,igsrc,ierfl,iptcmp)
  DIMENSION msrcs(numscs),ptsrec(2,numpts),irec(2,numpts),              &
            igsrc(numpts)
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  levelr = node(4,mrec)
!
  DO i=1,numpts
    igsrc(i)=0
  END DO
!
  DO lvlchk=levelr-1,levelr+1
!
    IF(lvlchk == 0) CYCLE

    DO ig=1,numscs
      msrc=msrcs(ig)
      IF(node(4,msrc) /= lvlchk) CYCLE
!
! now find points with msrc as source
!
      dx = rnode(9,msrc)
      dy = rnode(10,msrc)
      hx=2*hxposs(lvlchk)+xshift*dx
      hy=2*hyposs(lvlchk)+yshift*dy
      t1=rnode(16,msrc)-hx
      t2=rnode(17,msrc)+hx
      t3=rnode(18,msrc)-hy
      t4=rnode(19,msrc)+hy

      DO i=1,numpts
        ex1=ptsrec(1,i)*rnode(12,msrc)+ptsrec(2,i)*rnode(13,msrc)
        ex2=ptsrec(1,i)*rnode(14,msrc)+ptsrec(2,i)*rnode(15,msrc)
        IF((t1-ex1 >= 0).AND.(ex1-t2 >= 0).AND.                         &
            (t3-ex2 >= 0).AND.(ex2-t4 >= 0)) igsrc(i)=msrc
      END DO
    END DO

  END DO
!
! next check if source has been found for all points
!
  ierfl0=ierfl
  IF(ierfl == 1) THEN
    ierfl=0
    DO i=1,numpts
      IF(igsrc(i) == 0) ierfl=ierfl+1
    END DO

    IF(ierfl > 0) THEN
!
!   some points do not have a source, problems here
!
      WRITE(6,'(''  IN SOURCE, POINTS LACKING SOURCE, GRID '',i4)') mrec
      WRITE(6,'(''  POSSIBLE SOURCE GRIDS '',10I4)')(msrcs(i),i=1,numscs)
      DO i=1,numpts
        IF(igsrc(i) == 0) WRITE(6,240) ptsrec(1,i),ptsrec(2,i)
      END DO
240   FORMAT('  point  ',2(1X,e12.5))
      ierfl=1
    END IF

  END IF
!
!  finally, return only those points that have source
!
  IF(ierfl0 == 1.AND.ierfl == 0) RETURN
  IF(iptcmp == 1) THEN

    ist=0
    DO i=1,numpts
      IF(igsrc(i) /= 0)THEN
        ist=ist+1
        igsrc(ist)=igsrc(i)
        ptsrec(1,ist)=ptsrec(1,i)
        ptsrec(2,ist)=ptsrec(2,i)
        irec(1,ist)=irec(1,i)
        irec(2,ist)=irec(2,i)
      END IF
    END DO
    numpts=ist

    310   CONTINUE
  END IF

  RETURN
END SUBROUTINE getsrc

!
!--------------------------------------------------------------------------
!

SUBROUTINE calwht(mrec,msrc,points,isrc,whts,numpts,nwhts,              &
           ivar,ndim)
  DIMENSION points(2,numpts),whts(nwhts,numpts),isrc(2,numpts)
!
!  this routine calculates the weights for the interpolation
!  of mass points or vertical velocity (stuff that doesn't need
!  any rotation from grid msrc to grid mrec
!
!  if nw = 4 bilinear interp., = 9 bi-quadratic
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'

  fq(xp,x1,x2,x3,x4)=(xp-x2)*(xp-x3)/                                   &
                    ((x1-x2)*(x1-x3))+x4
  fl(xp,x1,x2)=(xp-x2)/(x1-x2)
!
  cosc = rnode( 21,msrc)
  sinc = rnode( 22,msrc)
  dxc  = rnode(9,msrc)
  dyc  = rnode(10,msrc)
  odxc = 1./dxc
  odyc = 1./dyc
  xor  = rnode(1,msrc)
  yor  = rnode(2,msrc)
!
  IF(ndim == 2) THEN
    xshift=stgxy(1,ivar)
    yshift=stgxy(2,ivar)
  ELSE
    xshift=stgxyz(1,ivar)
    yshift=stgxyz(2,ivar)
  END IF

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
  IF(nwhts == 9) THEN
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
END SUBROUTINE calwht
!
!--------------------------------------------------------------------------
!

SUBROUTINE hrzavg(vout,vin,nx,ny,nz)
!
! Perform horizontal plane average over nine points.
!
  DIMENSION vout(nx,ny,nz),vin(nx,ny,nz)

  a19 = 1.0/9.0
  DO k=1,nz
    DO j=2,ny-1
      DO i=2,nx-1
        vout(i,j,k)=(vin(i,j,k)+vin(i-1,j,k)+vin(i+1,j,k)               &
             +vin(i,j-1,k)+vin(i,j+1,k)+vin(i+1,j+1,k)                  &
             +vin(i+1,j-1,k)+vin(i-1,j+1,k)+vin(i-1,j-1,k))*a19
      END DO
    END DO
  END DO
  DO k=1,nz
    DO j=1,ny
      vout(1 ,j,k)=vin(1, j,k)
      vout(nx,j,k)=vin(nx,j,k)
    END DO
    DO i=1,nx
      vout(i,1 ,k)=vin(i,1 ,k)
      vout(i,ny,k)=vin(i,ny,k)
    END DO
  END DO
  RETURN
END SUBROUTINE hrzavg

  FUNCTION smstg2(ivar1,ivar2, ndim)
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'

  LOGICAL :: smstg2
  smstg2=.false.
  eps=1.0E-6
  CALL chkdim(ndim,'SMSTG2')
  IF(ndim == 2) THEN
    IF(ABS(stgxy (1,ivar1)-stgxy (1,ivar2 )) < eps .AND.                &
        ABS(stgxy (2,ivar1)-stgxy (2,ivar2)) < eps) smstg2=.true.
  ELSE
    IF(ABS(stgxyz(1,ivar1)-stgxyz(1,ivar2 )) < eps .AND.                &
        ABS(stgxyz(2,ivar1)-stgxyz(2,ivar2)) < eps) smstg2=.true.
  END IF
  RETURN
  END FUNCTION smstg2

  FUNCTION smstg3(ivar1,ivar2, ndim)
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  LOGICAL :: smstg3
  smstg3=.false.
  eps=1.0E-6
  CALL chkdim(ndim,'SMSTG3')
  IF(ndim == 2) THEN
    IF(ABS(stgxy (1,ivar1)-stgxy (1,ivar2)) < eps.AND.                  &
        ABS(stgxy (2,ivar1)-stgxy (2,ivar2)) < eps) smstg3=.true.
  ELSE
    IF(ABS(stgxyz(1,ivar1)-stgxyz(1,ivar2)) < eps.AND.                  &
        ABS(stgxyz(2,ivar1)-stgxyz(2,ivar2)) < eps.AND.                 &
        ABS(stgxyz(3,ivar1)-stgxyz(3,ivar2)) < eps) smstg3=.true.
  END IF
  RETURN
  END FUNCTION smstg3
!

SUBROUTINE addcst(mptr,ndim,ivar,cst)
!
! to add a constant CST to ndim-D variable IVAR for grid MPTR.
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'

  nx=node(5,mptr)
  ny=node(6,mptr)

  CALL chkdim(ndim,'GTRCRD')
  IF(ndim == 2) THEN
    nz=1
    iptr=igtnxy(mptr,ivar,1)
  ELSE
    nz=node(14,mptr)
    iptr=igtxyz(mptr,ivar,1)
  END IF

  DO ijk=1,nx*ny*nz
    a(ijk-1+iptr)=a(ijk-1+iptr)+cst
  END DO
  IF(ndim == 2) THEN
    CALL retnxy(mptr,ivar,1,iptr,.true.)
  ELSE
    CALL retxyz(mptr,ivar,1,iptr,.true.)
  END IF
  RETURN
END SUBROUTINE addcst
