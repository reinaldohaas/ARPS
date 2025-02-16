!
!--------------------------------------------------------------------------
!

SUBROUTINE updbc( mptr,samlvl,tmwght )
!
!  Provide boundary values (conditons) for grid MPTR by
!  interpolating from course grid or from overlapping fine grids
!  at the same level depending on the value of samlvl..
!
!  MPTR is the grid needing boundary values by interpolation
!  samlvl is a logical, if true  - interp from same level
!                              (i.e., from overlapping grid)
!                       false - interp from next level
!  TMWGHT is the weight for the interpolated value, with
!
!  phi(t)=TMWGHT*phi(interp)+(1-TMWGHT)*phi(t-dt(receive))   (1)
!
!  where phi(interp) is the interpolated value from source grid
!  at t+dt(source).
!
!  This routine resets phi(t-dt(receive)) at the boundaries so that
!  formula (1) is satisfied, then linear extrapolation of phi
!  from its value at t and t-dt(receive) will set phi(t+n*dt(receive))
!  equal to the interpolated value, where n=dt(source)/dt(receive)..
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricpu.inc'
  INCLUDE 'agricst.inc'
  DIMENSION msrcs(30),ibcvar(nxy2d+nxyz3d),ibcsta(nxy2d+nxyz3d)
  LOGICAL :: samlvl,smstg2
!
! here we assume the number of source grids never exceeds 30
!
  IF( intrpodr == 2 ) THEN
    nwhts=9
  ELSE
    nwhts=4
  END IF

  level = node(4,mptr)
  IF(level > lfine.OR.level == 0) RETURN
  cpu0 = f_cputime()
!
! Get the possible source grids for the b.c.'s of grid mptr.
!
  numscs=0
  lchk = level
  IF(.NOT.samlvl) lchk = lchk-1
!
!   do 20 lchk=max(1,level-1),level
!
  mgrid=lstart(lchk)
  15      IF(mgrid == 0) GO TO 20
  numscs=numscs+1
  msrcs(numscs)=mgrid
  mgrid=node(10,mgrid)
  GO TO 15
  20    CONTINUE
  5     CONTINUE

  nx=node(5,mptr)
  ny=node(6,mptr)
!
! --------------------------------------------------------------------
! calculate b.c.'s for 2-d scalars
! --------------------------------------------------------------------
!
  ndim=2
!
! initialize b.c. calculating status with value 1 for each scalar
!
  DO i=1,nxy2d
    ibcsta(i)=1
  END DO

  100   CONTINUE
!
! search for variables that need b.c.'s and group them according to
! their staggering
!
  nbcvar=0
  DO ibc=1,nxy2d
!
! look for 1st 2-d scalar that needs b.c. calculation
!
    IF(ibdxy(1,ibc) == 1.AND.ibcsta(ibc) == 1.AND. ibdxy(2,ibc) == 0)THEN
      nbcvar=nbcvar+1
      ibcvar(nbcvar)=ibc
      ibcsta(ibc)=0
      DO ivar=ibc+1,nxy2d
!
! list additonal variables with matching staggering to that of ibc.
!
        IF(ibdxy(1,ivar) == 1.AND.ibcsta(ivar) == 1.AND.                &
              ibdxy(2,ivar) == 0.AND.smstg2(ivar,ibc,2))THEN
          nbcvar=nbcvar+1
          ibcvar(nbcvar)=ivar
          ibcsta(ivar)=0
        END IF
      END DO
    END IF
  END DO
!
! if the b.c.'s for all 2-d scalars are calculated goto next step
!
  IF(nbcvar == 0)GO TO 200
!
! perform interpolations for 2-d scalar variables listed in ibcvar(nbcvar).
!
  mxpts=2*(nx+ny)
  iptr1=igetsp(2*mxpts)
  iptr2=igetsp(2*mxpts)
  iptr3=igetsp(2*mxpts)
  iptr4=igetsp(  mxpts)
  iptr5=igetsp(nwhts*mxpts)
  iptr6=igetsp(2*mxpts)
  iptr7=igetsp(2*mxpts)
!
  CALL itpbsl(mptr,msrcs,numscs,a(iptr1),a(iptr2),a(iptr3),             &
       a(iptr4),a(iptr5),a(iptr6),a(iptr7),                             &
       mxpts,nwhts,ibcvar,nbcvar,ndim,samlvl,tmwght)
!
!  restore tem space
  CALL resett
!
! go back to check if any more variable needs b.c. calculation
!
  GO TO 100
  200   CONTINUE
!
! --------------------------------------------------------------------
! calculate b.c.'s for 3-d scalars
! --------------------------------------------------------------------
!
  ndim=3
!
! initialize b.c. calculating status with value 1 for each vari.
!
  DO i=1,nxyz3d
    ibcsta(i)=1
  END DO
!
! search for variables that need b.c.'s and group them according to
! their staggering
!
  300   CONTINUE
  nbcvar=0
  DO ibc=1,nxyz3d
!
! look for 1st 3-d scalar that needs b.c. calculation
!
    IF(ibdxyz(1,ibc) == 1.AND.ibcsta(ibc) == 1.AND. ibdxyz(2,ibc) == 0)THEN
      nbcvar=nbcvar+1
      ibcvar(nbcvar)=ibc
      ibcsta(ibc)=0
      DO ivar=ibc+1,nxyz3d
!
! list additonal variables with matching staggering to that of ibc.
!
        IF(ibdxyz(1,ivar) == 1.AND.ibcsta(ivar) == 1.AND.               &
              ibdxyz(2,ivar) == 0.AND.smstg2(ivar,ibc,3))THEN
          nbcvar=nbcvar+1
          ibcvar(nbcvar)=ivar
          ibcsta(ivar)=0
        END IF
      END DO
    END IF
  END DO
!
! if the b.c.'s for all 3-d scalars are calculated goto next step
!
  IF(nbcvar == 0) GO TO 400
!
! perform interpolations for 3-d scalars listed in ibcvar(nbcvar).
!
  mxpts=2*(nx+ny)
  iptr1=igetsp(2*mxpts)
  iptr2=igetsp(2*mxpts)
  iptr3=igetsp(2*mxpts)
  iptr4=igetsp(  mxpts)
  iptr5=igetsp(nwhts*mxpts)
  iptr6=igetsp(2*mxpts)
  iptr7=igetsp(2*mxpts)
!
  CALL itpbsl(mptr,msrcs,numscs,a(iptr1),a(iptr2),a(iptr3),             &
       a(iptr4),a(iptr5),a(iptr6),a(iptr7),                             &
       mxpts,nwhts,ibcvar,nbcvar,ndim,samlvl,tmwght)
!
!  restore tem space
!
  CALL resett
!
! go back to check if any more 3-d scalars needs b.c. calculation
!
  GO TO 300
  400   CONTINUE
!
!-------------------------------------------------------------------
! calculate b.c.'s for 2-d vectors
! --------------------------------------------------------------------
!
  ndim=2
!
! initialize b.c. calculating status with value 1 for each vari.
!
  DO i=1,nxy2d
    ibcsta(i)=1
  END DO

  500   CONTINUE
!
! search for 2-d vectors that need b.c.'s
!
  nbcvar=0
  DO ibc=1,nxy2d
!
! look for 1st 2-d vector that needs b.c. calculation
!
    IF(ibdxy(1,ibc) == 1.AND.ibdxy(2,ibc) > 0.AND. ibcsta(ibc) == 1) THEN
      ibc2=ibdxy(2,ibc)
      ibcvar(nbcvar+1)=ibc
      ibcvar(nbcvar+2)=ibc2
      nbcvar=nbcvar+2
      ibcsta(ibc )=0
      ibcsta(ibc2)=0
      EXIT
    END IF
  END DO
  67    CONTINUE
!
! if the b.c.'s for all 2-d vectors are calculated goto next step
!
  IF(nbcvar == 0)GO TO 600
!
! perform interpolations for 2-d vectors listed in ibcvar(nbcvar).
!
  mxpts=2*(nx+ny)
  iptr1=igetsp(4*mxpts)
  iptr2=igetsp(4*mxpts)
  iptr3=igetsp(4*mxpts)
  iptr4=igetsp(2*mxpts)
  iptr5=igetsp(2*nwhts*mxpts)
  iptr6=igetsp(4*mxpts)
  iptr7=igetsp(4*mxpts)
!
  CALL itpbvt(mptr,msrcs,numscs,a(iptr1),a(iptr2),a(iptr3),             &
       a(iptr4),a(iptr5),a(iptr6),a(iptr7),                             &
       mxpts,nwhts,ibcvar,nbcvar,ndim,samlvl,tmwght)
!
!  restore tem space
!
  CALL resett
!
! go back to check if any more 2-d vector needs b.c. calculation
!
  GO TO 500
  600   CONTINUE
!
!------------------------------------------------------------------
! calculate b.c.'s for 3-d vectors
!------------------------------------------------------------------
!
  ndim=3
!
! initialize b.c. calculating status with value 1 for each vari.
!
  DO i=1,nxyz3d
    ibcsta(i)=1
  END DO
!
  700   CONTINUE
!
! search for variables that need b.c.'s and group them according to
! their staggering
!
  nbcvar=0
  DO ibc=1,nxyz3d
!
! look for 1st 3-d vector that needs b.c. calculation
!
    IF(ibdxyz(1,ibc) == 1.AND.ibdxyz(2,ibc) > 0.AND. ibcsta(ibc) == 1) THEN
      ibc2=ibdxyz(2,ibc)
      ibcvar(nbcvar+1)=ibc
      ibcvar(nbcvar+2)=ibc2
      nbcvar=nbcvar+2
      ibcsta(ibc )=0
      ibcsta(ibc2)=0
      EXIT
    END IF
  END DO
  92    CONTINUE
!
! if the b.c.'s for all 3-d vectors are calculated goto next step
!
  IF(nbcvar == 0)GO TO 800
!
! perform interpolations for 3-d vectors listed in ibcvar(nbcvar).
!
  mxpts=2*(nx+ny)
  iptr1=igetsp(4*mxpts)
  iptr2=igetsp(4*mxpts)
  iptr3=igetsp(4*mxpts)
  iptr4=igetsp(2*mxpts)
  iptr5=igetsp(2*nwhts*mxpts)
  iptr6=igetsp(4*mxpts)
  iptr7=igetsp(4*mxpts)
!
  IF( verbose6 ) WRITE(6,'('' ITPBVT CALLED FOR 3-D VECTORS'')')
!
  CALL itpbvt(mptr,msrcs,numscs,a(iptr1),a(iptr2),a(iptr3),             &
       a(iptr4),a(iptr5),a(iptr6),a(iptr7),                             &
       mxpts,nwhts,ibcvar,nbcvar,ndim,samlvl,tmwght)
!
  IF( verbose6 ) THEN
    DO ib=1,nbcvar,2
      WRITE(6,'(2i5)') ibcvar(ib),ibcvar(ib+1)
    END DO
  END IF
!
!  restore tem space
!
  CALL resett
!
! go back to check if any more 3-d vector needs b.c. calculation
!
  GO TO 700
  800   CONTINUE
  cpu_bndcint = cpu_bndcint + f_cputime() - cpu0
  RETURN
END SUBROUTINE updbc
!
!--------------------------------------------------------------------------
!

SUBROUTINE itpbsl(mrec,msrcs,numscs,ptsrec,                             &
           irec,isrc,igsrc,whts,pts,irecp,mxpts,nwhts,                  &
           ibcvar,nbcvar,ndim,samlvl,tmwght)
!
  DIMENSION msrcs(numscs),ptsrec(2,mxpts),irec(2,mxpts),                &
            isrc(2,mxpts),igsrc(mxpts),whts(nwhts,mxpts),               &
            pts(2,mxpts),irecp(2,mxpts),                                &
            ibcvar(nbcvar)

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  LOGICAL :: samlvl
  nw1=nint( SQRT(FLOAT(nwhts)) )
  nw2=nw1
!
! Calculate and store the coor's of p-points of the receive
! grid in array ptsrec with their indecies stored in irec.
!
  IF(nbcvar == 0) RETURN
  CALL chkdim(ndim,'ITPBSL')

  nx=node(5,mrec)
  ny=node(6,mrec)
  cosr=rnode(21,mrec)
  sinr=rnode(22,mrec)
  dxr=rnode(9 ,mrec)
  dyr=rnode(10,mrec)
  xor=rnode(1,mrec)
  yor=rnode(2,mrec)

  ivar=ibcvar(1)
  IF(ndim == 2) THEN
    xshift=stgxy(1,ivar)
    yshift=stgxy(2,ivar)
    nx1=nx-idmxy(1,ivar)
    ny1=ny-idmxy(2,ivar)
  ELSE
    xshift=stgxyz(1,ivar)
    yshift=stgxyz(2,ivar)
    nx1=nx-idmxyz(1,ivar)
    ny1=ny-idmxyz(2,ivar)
  END IF

  xc=xshift*dxr*cosr
  ys=yshift*dyr*sinr
  xs=xshift*dxr*sinr
  yc=yshift*dyr*cosr

  xorp=xor+xc-ys
  yorp=yor+xs+yc
!
! To calculate the absolute coordinates of the model boundary points
!
! Attention:
!
! Here we assume that in x direction the boundaries are at i=1 and nx1.
! and at j=1 and ny1 in y direction.
! If in a particular solver (model), the boundaries for certain
! variables are defined at points other the above, the use needs to
! change the following loop accordingly (i.e. change ist, jst, iend
! and jend). As a result, the part in routine UPDSCL that calculates
! the points for interior points updating also need change.
!
  ist=1
  jst=1
  iend=nx1
  jend=ny1
  DO j=jst,jend,jend-jst
    jj=(j-jst)/(jend-jst)*(iend-ist+1)
    DO i=ist,iend
      ii=jj+i-ist+1
      x=(i-1)*dxr
      y=(j-1)*dyr
      ptsrec(1,ii) = xorp+x*cosr-y*sinr
      ptsrec(2,ii) = yorp+x*sinr+y*cosr
      irec  (1,ii) = i
      irec  (2,ii) = j
    END DO
  END DO
  numpts=MAX(0, 2*(iend-ist+1))

  DO i=ist,iend,iend-ist
    ii=(i-ist)/(iend-ist)*(jend-jst-1)
    DO j=jst+1,jend-1
      jj=ii+j-jst+numpts
      x=(i-1)*dxr
      y=(j-1)*dyr
      ptsrec(1,jj) = xorp+x*cosr-y*sinr
      ptsrec(2,jj) = yorp+x*sinr+y*cosr
      irec  (1,jj) = i
      irec  (2,jj) = j
    END DO
  END DO
  numpts=numpts+MAX(0,2*(jend-jst-1))
  npts0=numpts
!
! Find the source points for receive grid points ptsrec
!
  IF(samlvl) THEN
    ierfl=0
    iptcmp=1
  ELSE
    ierfl=1
    iptcmp=1
  END IF
  CALL getsrc(mrec,msrcs,numscs,xshift,yshift,ptsrec,numpts,            &
       irec,igsrc,ierfl,iptcmp)
!
! if no point can draw boudary value from listed source grids, skip through
!
  IF(.NOT.samlvl) THEN
    IF(numpts /= npts0) THEN
      WRITE(6,'(i4,'' BOUNDARY POINTS DID NOT FIND SOURCE FOR B.C.S''   &
      &   )') npts0-numpts
      WRITE(6,'(''IT WAS FOR '',I1,''-D SCALAR NO. '',I4)')ndim,ivar
      WRITE(6,'(''JOB STOPPED IN ITPBSL.'')')
      STOP
    END IF
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
!
! Interpolate for the boundary values for variables in list ibcvar
!
    DO i=1,nbcvar
      ivar=ibcvar(i)
      IF(ndim == 2) THEN
        n3rd=1
        irptr=igtnxy(mrec,ivar,1)
        istmdt=igtnxy(mrec,ibdxy(3,ivar),1)
        isptr=igtnxy(msrc,ivar,1)
      ELSE
        n3rd=nzs
        irptr=igtxyz(mrec,ivar,1)
        istmdt=igtxyz(mrec,ibdxyz(3,ivar),1)
        isptr=igtxyz(msrc,ivar,1)
      END IF

      CALL intpst1(a(irptr),a(isptr),a(istmdt),tmwght,                  &
          nxr,nyr,nxs,nys,n3rd,                                         &
          irecp,isrc,whts,npts,nw1,nw2)

      IF(ndim == 2) THEN
        CALL retnxy(mrec,ivar,1,irptr,.true.)
        CALL retnxy(mrec,ibdxy(3,ivar),1,istmdt,.false.)
        CALL retnxy(msrc,ivar,1,isptr,.false.)
      ELSE
        CALL retxyz(mrec,ivar,1,irptr,.false.)
        CALL retxyz(mrec,ibdxyz(3,ivar),1,istmdt,.true.)
!      call retxyz(mrec,ivar,1,irptr,.true.)
!      call retxyz(mrec,ibdxyz(3,ivar),1,istmdt,.false.)
        CALL retxyz(msrc,ivar,1,isptr,.false.)
      END IF
    END DO
!
  END DO
  500   CONTINUE
  RETURN
END SUBROUTINE itpbsl
!
!--------------------------------------------------------------------------
!

SUBROUTINE itpbvt(mrec,msrcs,numscs,ptsrec,                             &
           irec,isrc,igsrc,whts,pts,irecp,mxpts,nwhts,                  &
           ibcvar,nbcvar,ndim,samlvl,tmwght)
!
  DIMENSION msrcs(numscs),ptsrec(2,mxpts,2),irec(2,mxpts,2),            &
            isrc(2,mxpts,2),igsrc(mxpts,2),whts(nwhts,mxpts,2),         &
            pts(2,mxpts,2),irecp(2,mxpts,2),ibcvar(nbcvar),             &
            numpts(2),npts(2)

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  LOGICAL :: samstg, smstg2, samlvl

  nw1=nint( SQRT(FLOAT(nwhts)) )
  nw2=nw1
!
  IF(nbcvar == 0) RETURN
  CALL chkdim(ndim,'ITPBVT')

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
! To loop through all vectors, which has two components stored consecutively
! in nbcvar. The 1st compt. has to be the one in x direction.
!
  DO ivtr=1,nbcvar,2

    ivar1=ibcvar(ivtr)
    ivar2=ibcvar(ivtr+1)

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
! them.If the two compnts have the same staggering, they have the
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

      ivar=ibcvar(ivtr+icmpr-1)

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
! To calculate the absolute coordinates of the model boundary points
!
! Attention:
!
! Here we assume that in x direction the boundaries are at i=1 and nx1.
! and at j=1 and ny1 in y direction.
! If in a particular solver (model), the boundaries for certain
! variables are defined at points other the above, the use needs to
! change the following loop accordingly (i.e. change ist, jst, iend
! and jend). As a result, the part in routine UPDVTR that calculates
! the points for interior points updating also need change.
!
      ist=1
      jst=1
      iend=nx1
      jend=ny1
      DO j=jst,jend,jend-jst
        jj=(j-jst)/(jend-jst)*(iend-ist+1)
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
      numpts(icmpr)=MAX(0,2*(iend-ist+1))

      DO i=ist,iend,iend-ist
        ii=(i-ist)/(iend-ist)*(jend-jst-1)
        DO j=jst+1,jend-1
          jj=ii+j-jst+numpts(icmpr)
          x=(i-1)*dxr
          y=(j-1)*dyr
          ptsrec(1,jj,icmpr) = xorp+x*cosr-y*sinr
          ptsrec(2,jj,icmpr) = yorp+x*sinr+y*cosr
          irec  (1,jj,icmpr) = i
          irec  (2,jj,icmpr) = j
        END DO
      END DO
      numpts(icmpr)=numpts(icmpr)+MAX(0,2*(jend-jst-1))
      npts0=numpts(icmpr)
!
! Find the source grid for receive grid points ptsrec
!
      IF(samlvl) THEN
        ierfl=0
        iptcmp=1
      ELSE
        ierfl=1
        iptcmp=1
      END IF
      CALL getsrc(mrec,msrcs,numscs,xshft,yshft,ptsrec(1,1,icmpr),      &
          numpts(icmpr),irec(1,1,icmpr),igsrc(1,icmpr),ierfl,iptcmp)
!
      IF(.NOT.samlvl) THEN
        IF(numpts(icmpr) /= npts0) THEN
          WRITE(6,'(i4,'' BOUNDARY POINTS DID NOT FIND SOURCE FOR B.C.S'' &
          &   )') npts0-numpts(icmpr)
          WRITE(6,'(''IT WAS FOR '',I1,''-D VECTOR COMPNT. '',2I4)')    &
              ndim,ivar
          WRITE(6,'(''JOB STOPPED IN ITPBVT.'')')
          STOP
        END IF
      END IF
    END DO
!
! find storage pointers for the receive grid variables (2 compnts of vector).
! uppack them if neccesary.
!
    IF(ndim == 2) THEN
      irptr1=igtnxy(mrec,ivar1,1)
      irptr2=igtnxy(mrec,ivar2,1)
      mrptr1=igtnxy(mrec,ibdxy(3,ivar1),1)
      mrptr2=igtnxy(mrec,ibdxy(3,ivar2),1)
    ELSE
      irptr1=igtxyz(mrec,ivar1,1)
      irptr2=igtxyz(mrec,ivar2,1)
      mrptr1=igtxyz(mrec,ibdxyz(3,ivar1),1)
      mrptr2=igtxyz(mrec,ibdxyz(3,ivar2),1)
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
          ivsrc=ibcvar(ivtr+icmps-1)
          CALL calwht(mrec,msrc,pts(1,1,icmpr),isrc(1,1,icmpr),         &
                      whts(1,1,icmpr),npts(icmpr),nwhts,ivsrc,ndim)
        END DO
        IF(npts(1) == 0.AND.npts(2) == 0) CYCLE
!
! find pointer for one of the components (icmps) of the source grid
! vector. If it is packed, do uppacking for this variable.
!
        ivsrc=ibcvar(ivtr+icmps-1)
        IF(ndim == 2) THEN
          n3rd=1
          isptr=igtnxy(msrc,ivsrc,1)
        ELSE
          n3rd=nzs
          isptr=igtxyz(msrc,ivsrc,1)
        END IF

        IF(icmps == 1)THEN
          proj1=cossr
          proj2=-sinsr
        ELSE
          proj1=sinsr
          proj2=cossr
        END IF

        icmpr=1
        IF(npts(icmpr) /= 0) THEN
          CALL intpvt1(a(irptr1),a(isptr),a(mrptr1),tmwght,             &
              nxr,nyr,nxs,nys,n3rd,                                     &
              irecp(1,1,icmpr),isrc(1,1,icmpr),whts(1,1,icmpr),         &
              npts(icmpr),nw1,nw2,icmps,proj1)
        END IF

        icmpr=2
        IF(npts(icmpr) /= 0) THEN
          CALL intpvt1(a(irptr2),a(isptr),a(mrptr2),tmwght,             &
              nxr,nyr,nxs,nys,n3rd,                                     &
              irecp(1,1,icmpr),isrc(1,1,icmpr),whts(1,1,icmpr),         &
              npts(icmpr),nw1,nw2,icmps,proj2)
        END IF

        IF(ndim == 2) THEN
          CALL retnxy(msrc,ivsrc,1,isptr,.false.)
        ELSE
          CALL retxyz(msrc,ivsrc,1,isptr,.false.)
        END IF
      END DO
!
    END DO
    IF(ndim == 2) THEN
      CALL retnxy(mrec,ivar1,1,irptr1,.false.)
      CALL retnxy(mrec,ivar2,1,irptr2,.false.)
      CALL retnxy(mrec,ibdxy(3,ivar1),1,mrptr1,.true.)
      CALL retnxy(mrec,ibdxy(3,ivar2),1,mrptr2,.true.)
    ELSE
      CALL retxyz(mrec,ivar1,1,irptr1,.false.)
      CALL retxyz(mrec,ivar2,1,irptr2,.false.)
      CALL retxyz(mrec,ibdxyz(3,ivar1),1,mrptr1,.true.)
      CALL retxyz(mrec,ibdxyz(3,ivar2),1,mrptr2,.true.)
    END IF
  END DO
  RETURN
END SUBROUTINE itpbvt
!
!--------------------------------------------------------------------------
!

SUBROUTINE intpst(varrec,varsrc,varmdt,tmwght,                          &
           nxr,nyr,nxs,nys,nz,                                          &
           irec,isrc,whts,numpts,nw1,nw2)
  DIMENSION varrec(nxr,nyr,nz),varsrc(nxs,nys,nz),                      &
            varmdt(nxr,nyr,nz),                                         &
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
    1114    FORMAT(1X,4I7,4E10.3)
  END IF
!
!  set weights for temporal linear interp
!
  w1 = tmwght
  w2 = 1.-w1
!
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
          varrec(ir,jr,ks)=w1*(varsrc(is,js,ks)*whts(1,1,np)            &
              +varsrc(is+1,js,ks)*whts(2,1,np)                          &
              +varsrc(is,js+1,ks)*whts(1,2,np)                          &
              +varsrc(is+1,js+1,ks)*whts(2,2,np))                       &
              +w2*varmdt(ir,jr,ks)
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
          varrec(ir,jr,ks)=w1*(varsrc(is,js,ks)*whts(1,1,np)            &
              +varsrc(is+1,js,ks)*whts(2,1,np)                          &
              +varsrc(is,js+1,ks)*whts(1,2,np)                          &
              +varsrc(is+1,js+1,ks)*whts(2,2,np))                       &
              +w2*varmdt(ir,jr,ks)
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
          varrec(ir,jr,ks)=   w1*(                                      &
              varsrc(is-1,js-1,ks)*whts(1,1,np)                         &
              +varsrc(is  ,js-1,ks)*whts(2,1,np)                        &
              +varsrc(is+1,js-1,ks)*whts(3,1,np)                        &
              +varsrc(is-1,js  ,ks)*whts(1,2,np)                        &
              +varsrc(is  ,js  ,ks)*whts(2,2,np)                        &
              +varsrc(is+1,js  ,ks)*whts(3,2,np)                        &
              +varsrc(is-1,js+1,ks)*whts(1,3,np)                        &
              +varsrc(is  ,js+1,ks)*whts(2,3,np)                        &
              +varsrc(is+1,js+1,ks)*whts(3,3,np)  )                     &
              +w2*varmdt(ir,jr,ks)
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
          varrec(ir,jr,ks)=   w1*(                                      &
              varsrc(is-1,js-1,ks)*whts(1,1,np)                         &
              +varsrc(is  ,js-1,ks)*whts(2,1,np)                        &
              +varsrc(is+1,js-1,ks)*whts(3,1,np)                        &
              +varsrc(is-1,js  ,ks)*whts(1,2,np)                        &
              +varsrc(is  ,js  ,ks)*whts(2,2,np)                        &
              +varsrc(is+1,js  ,ks)*whts(3,2,np)                        &
              +varsrc(is-1,js+1,ks)*whts(1,3,np)                        &
              +varsrc(is  ,js+1,ks)*whts(2,3,np)                        &
              +varsrc(is+1,js+1,ks)*whts(3,3,np)  )                     &
              +w2*varmdt(ir,jr,ks)
        END DO
      END DO
    END IF

  END IF

  RETURN
END SUBROUTINE intpst
!
!--------------------------------------------------------------------------
!

SUBROUTINE intpvt(varrec,varsrc,varmdt,tmwght,                          &
           nxr,nyr,nxs,nys,nz,                                          &
           irec,isrc,whts,numpts,nw1,nw2,mode,proj)
  DIMENSION varrec(nxr,nyr,nz),varsrc(nxs,nys,nz),                      &
            varmdt(nxr,nyr,nz),                                         &
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
    1114    FORMAT(1X,4I7,4E10.3)
  END IF

  iaccum=0
  IF(mode == 2) iaccum=1
  w1 = tmwght
  w2 = 1.-w1

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
              +varsrc(is+1,js+1,ks)*whts(2,2,np))*proj*w1               &
              +w2*iaccum*varmdt(ir,jr,ks)
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
              +varsrc(is+1,js+1,ks)*whts(2,2,np))*proj*w1               &
              +w2*iaccum*varmdt(ir,jr,ks)
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
              +varsrc(is+1,js+1,ks)*whts(3,3,np))*proj*w1               &
              +w2*iaccum*varmdt(ir,jr,ks)
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
              +varsrc(is+1,js+1,ks)*whts(3,3,np))*proj*w1               &
              +w2*iaccum*varmdt(ir,jr,ks)
        END DO
      END DO
    END IF

  END IF

  RETURN
END SUBROUTINE intpvt
!
!--------------------------------------------------------------------------
!

SUBROUTINE intpst1(varrec,varsrc,varmdt,tmwght,                         &
           nxr,nyr,nxs,nys,nz,                                          &
           irec,isrc,whts,numpts,nw1,nw2)
  DIMENSION varrec(nxr,nyr,nz),varsrc(nxs,nys,nz),                      &
            varmdt(nxr,nyr,nz),                                         &
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
    1114    FORMAT(1X,4I7,4E10.3)
  END IF
!
!  set weights for temporal linear interp
!
  w1 = tmwght
  w2 = 1.-w1
!
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
          varmdt(ir,jr,ks)=(varrec(ir,jr,ks)-w1*                        &
              (varsrc(is,js,ks)*whts(1,1,np)                            &
              +varsrc(is+1,js,ks)*whts(2,1,np)                          &
              +varsrc(is,js+1,ks)*whts(1,2,np)                          &
              +varsrc(is+1,js+1,ks)*whts(2,2,np)))/w2
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
          varmdt(ir,jr,ks)=(varrec(ir,jr,ks)-w1*                        &
              (varsrc(is,js,ks)*whts(1,1,np)                            &
              +varsrc(is+1,js,ks)*whts(2,1,np)                          &
              +varsrc(is,js+1,ks)*whts(1,2,np)                          &
              +varsrc(is+1,js+1,ks)*whts(2,2,np)))/w2
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
          varmdt(ir,jr,ks)=(varrec(ir,jr,ks)-w1*(                       &
              varsrc(is-1,js-1,ks)*whts(1,1,np)                         &
              +varsrc(is  ,js-1,ks)*whts(2,1,np)                        &
              +varsrc(is+1,js-1,ks)*whts(3,1,np)                        &
              +varsrc(is-1,js  ,ks)*whts(1,2,np)                        &
              +varsrc(is  ,js  ,ks)*whts(2,2,np)                        &
              +varsrc(is+1,js  ,ks)*whts(3,2,np)                        &
              +varsrc(is-1,js+1,ks)*whts(1,3,np)                        &
              +varsrc(is  ,js+1,ks)*whts(2,3,np)                        &
              +varsrc(is+1,js+1,ks)*whts(3,3,np)))/w2
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
          varmdt(ir,jr,ks)=(varrec(ir,jr,ks)-w1*(                       &
              varsrc(is-1,js-1,ks)*whts(1,1,np)                         &
              +varsrc(is  ,js-1,ks)*whts(2,1,np)                        &
              +varsrc(is+1,js-1,ks)*whts(3,1,np)                        &
              +varsrc(is-1,js  ,ks)*whts(1,2,np)                        &
              +varsrc(is  ,js  ,ks)*whts(2,2,np)                        &
              +varsrc(is+1,js  ,ks)*whts(3,2,np)                        &
              +varsrc(is-1,js+1,ks)*whts(1,3,np)                        &
              +varsrc(is  ,js+1,ks)*whts(2,3,np)                        &
              +varsrc(is+1,js+1,ks)*whts(3,3,np)))/w2
        END DO
      END DO
    END IF

  END IF

  RETURN
END SUBROUTINE intpst1
!
!--------------------------------------------------------------------------
!

SUBROUTINE intpvt1(varrec,varsrc,varmdt,tmwght,                         &
           nxr,nyr,nxs,nys,nz,                                          &
           irec,isrc,whts,numpts,nw1,nw2,mode,proj)
  DIMENSION varrec(nxr,nyr,nz),varsrc(nxs,nys,nz),                      &
            varmdt(nxr,nyr,nz),                                         &
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
    1114    FORMAT(1X,4I7,4E10.3)
  END IF

  iaccum=0
  IF(mode == 2) iaccum=1
  w1 = tmwght
  w2 = 1.-w1

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
          varmdt(ir,jr,ks)=varmdt(ir,jr,ks)*iaccum                      &
              +(varrec(ir,jr,ks)*iaccum                                 &
              -(varsrc(is,js,ks)*whts(1,1,np)                           &
              +varsrc(is+1,js,ks)*whts(2,1,np)                          &
              +varsrc(is,js+1,ks)*whts(1,2,np)                          &
              +varsrc(is+1,js+1,ks)*whts(2,2,np))*proj*w1)/w2
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
          varmdt(ir,jr,ks)=varmdt(ir,jr,ks)*iaccum                      &
              +(varrec(ir,jr,ks)*iaccum                                 &
              -(varsrc(is,js,ks)*whts(1,1,np)                           &
              +varsrc(is+1,js,ks)*whts(2,1,np)                          &
              +varsrc(is,js+1,ks)*whts(1,2,np)                          &
              +varsrc(is+1,js+1,ks)*whts(2,2,np))*proj*w1)/w2
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
          varmdt(ir,jr,ks)=varmdt(ir,jr,ks)*iaccum                      &
              +(varrec(ir,jr,ks)*iaccum                                 &
              -(varsrc(is-1,js-1,ks)*whts(1,1,np)                       &
              +varsrc(is  ,js-1,ks)*whts(2,1,np)                        &
              +varsrc(is+1,js-1,ks)*whts(3,1,np)                        &
              +varsrc(is-1,js  ,ks)*whts(1,2,np)                        &
              +varsrc(is  ,js  ,ks)*whts(2,2,np)                        &
              +varsrc(is+1,js  ,ks)*whts(3,2,np)                        &
              +varsrc(is-1,js+1,ks)*whts(1,3,np)                        &
              +varsrc(is  ,js+1,ks)*whts(2,3,np)                        &
              +varsrc(is+1,js+1,ks)*whts(3,3,np))*proj*w1)/w2
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
          varmdt(ir,jr,ks)=varmdt(ir,jr,ks)*iaccum                      &
              +(varrec(ir,jr,ks)*iaccum                                 &
              -(varsrc(is-1,js-1,ks)*whts(1,1,np)                       &
              +varsrc(is  ,js-1,ks)*whts(2,1,np)                        &
              +varsrc(is+1,js-1,ks)*whts(3,1,np)                        &
              +varsrc(is-1,js  ,ks)*whts(1,2,np)                        &
              +varsrc(is  ,js  ,ks)*whts(2,2,np)                        &
              +varsrc(is+1,js  ,ks)*whts(3,2,np)                        &
              +varsrc(is-1,js+1,ks)*whts(1,3,np)                        &
              +varsrc(is  ,js+1,ks)*whts(2,3,np)                        &
              +varsrc(is+1,js+1,ks)*whts(3,3,np))*proj*w1)/w2
        END DO
      END DO
    END IF

  END IF

  RETURN
END SUBROUTINE intpvt1
!
!--------------------------------------------------------------------------
!

SUBROUTINE exchbc( mptr1,mptr2 )
!
!  Exchange boundary values (conditons) by interpolation bwtween
!  two overlapping grids at the same level. mptr1 and mptr2 are the
!  grid number.
!
!  Calling this once instead of calling updbc(mptr) twice so that each
!  grid do not need to be unpacked and packed twice.
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
  DIMENSION ibcvar(nxy2d+nxyz3d),ibcsta(nxy2d+nxyz3d)
  LOGICAL :: smstg2
!
  IF( intrpodr == 2 ) THEN
    nwhts=9
  ELSE
    nwhts=4
  END IF

  level1 = node(4,mptr1)
  IF(level1 > lfine.OR.level1 == 0) GO TO 999
  level2 = node(4,mptr2)
  IF(level2 > lfine.OR.level2 == 0) GO TO 998
  IF(level1 /= level2) GO TO 997
  nx1=node(5,mptr1)
  ny1=node(6,mptr1)
  nx2=node(5,mptr2)
  ny2=node(6,mptr2)
!
! --------------------------------------------------------------------
! calculate b.c.'s for 2-d scalars
! --------------------------------------------------------------------
!
  ndim=2
!
! initialize b.c. calculating status with value 1 for each scalar
!
  DO i=1,nxy2d
    ibcsta(i)=1
  END DO

  100   CONTINUE
!
! search for variables that need b.c.'s and group them according to
! their staggering
!
  nbcvar=0
  DO ibc=1,nxy2d
!
! look for 1st 2-d scalar that needs b.c. calculation
!
    IF(ibdxy(1,ibc) == 1.AND.ibcsta(ibc) == 1.AND. ibdxy(2,ibc) == 0)THEN
      nbcvar=nbcvar+1
      ibcvar(nbcvar)=ibc
      ibcsta(ibc)=0
      DO ivar=ibc+1,nxy2d
!
! list additonal variables with matching staggering to that of ibc.
!
        IF(ibdxy(1,ivar) == 1.AND.ibcsta(ivar) == 1.AND.                &
              ibdxy(2,ivar) == 0.AND.smstg2(ivar,ibc,ndim))THEN
          nbcvar=nbcvar+1
          ibcvar(nbcvar)=ivar
          ibcsta(ivar)=0
        END IF
      END DO
    END IF
  END DO
!
! if the b.c.'s for all 2-d scalars are calculated goto next step
!
  IF(nbcvar == 0)GO TO 200
!
! perform interpolations for 2-d scalar variables listed in ibcvar(nbcvar).
!
  mxpts=2*MAX( nx1+ny1,nx2+ny2 )
  iptr1=igetsp(2*mxpts)
  iptr2=igetsp(2*mxpts*2)
  iptr3=igetsp(2*mxpts*2)
  iptr4=igetsp(  mxpts)
  iptr5=igetsp(nwhts*mxpts*2)
  iptr6=igetsp(2*mxpts)
!
  CALL itpbs2(mptr1,mptr2,a(iptr1),a(iptr2),a(iptr3),a(iptr4),          &
       a(iptr5),a(iptr6),mxpts,nwhts,ibcvar,nbcvar,ndim)
!
!  restore tem space
  CALL resett
!
! go back to check if any more variable needs b.c. calculation
!
  GO TO 100
  200   CONTINUE
!
! --------------------------------------------------------------------
! calculate b.c.'s for 3-d scalars
! --------------------------------------------------------------------
!
  ndim=3
!
! initialize b.c. calculating status with value 1 for each vari.
!
  DO i=1,nxyz3d
    ibcsta(i)=1
  END DO
!
! search for variables that need b.c.'s and group them according to
! their staggering
!
  300   CONTINUE
  nbcvar=0
  DO ibc=1,nxyz3d
!
! look for 1st 3-d scalar that needs b.c. calculation
!
    IF(ibdxyz(1,ibc) == 1.AND.ibcsta(ibc) == 1.AND. ibdxyz(2,ibc) == 0)THEN
      nbcvar=nbcvar+1
      ibcvar(nbcvar)=ibc
      ibcsta(ibc)=0
      DO ivar=ibc+1,nxyz3d
!
! list additonal variables with matching staggering to that of ibc.
!
        IF(ibdxyz(1,ivar) == 1.AND.ibcsta(ivar) == 1.AND.               &
              ibdxyz(2,ivar) == 0.AND.smstg2(ivar,ibc,ndim))THEN
          nbcvar=nbcvar+1
          ibcvar(nbcvar)=ivar
          ibcsta(ivar)=0
        END IF
      END DO
    END IF
  END DO
!
! if the b.c.'s for all 3-d scalars are calculated goto next step
!
  IF(nbcvar == 0) GO TO 400
!
! perform interpolations for 3-d scalars listed in ibcvar(nbcvar).
!
  mxpts=2*MAX( nx1+ny1,nx2+ny2 )
  iptr1=igetsp(2*mxpts)
  iptr2=igetsp(2*mxpts*2)
  iptr3=igetsp(2*mxpts*2)
  iptr4=igetsp(  mxpts)
  iptr5=igetsp(nwhts*mxpts*2)
  iptr6=igetsp(2*mxpts)
!
  CALL itpbs2(mptr1,mptr2,a(iptr1),a(iptr2),a(iptr3),a(iptr4),          &
       a(iptr5),a(iptr6),mxpts,nwhts,ibcvar,nbcvar,ndim)
!
!  restore tem space
!
  CALL resett
!
! go back to check if any more 3-d scalars needs b.c. calculation
!
  GO TO 300
  400   CONTINUE
!
!-------------------------------------------------------------------
! calculate b.c.'s for 2-d vectors
! --------------------------------------------------------------------
!
  ndim=2
!
! initialize b.c. calculating status with value 1 for each vari.
!
  DO i=1,nxy2d
    ibcsta(i)=1
  END DO

  500   CONTINUE
!
! search for 2-d vectors that need b.c.'s
!
  nbcvar=0
  DO ibc=1,nxy2d
!
! look for 1st 2-d vector that needs b.c. calculation
!
    IF(ibdxy(1,ibc) == 1.AND.ibdxy(2,ibc) > 0.AND. ibcsta(ibc) == 1) THEN
      ibc2=ibdxy(2,ibc)
      ibcvar(nbcvar+1)=ibc
      ibcvar(nbcvar+2)=ibc2
      nbcvar=nbcvar+2
      ibcsta(ibc )=0
      ibcsta(ibc2)=0
      EXIT
    END IF
  END DO
  67    CONTINUE
!
! if the b.c.'s for all 2-d vectors are calculated goto next step
!
  IF(nbcvar == 0)GO TO 600
!
! perform interpolations for 2-d vectors listed in ibcvar(nbcvar).
!
  mxpts=2*MAX( nx1+ny1,nx2+ny2 )
  iptr1=igetsp(4*mxpts)
  iptr2=igetsp(4*mxpts)
  iptr3=igetsp(4*mxpts)
  iptr4=igetsp(2*mxpts)
  iptr5=igetsp(nwhts*mxpts*2)
  iptr6=igetsp(4*mxpts)
!
  CALL itpbv2(mptr1,mptr2,a(iptr1),a(iptr2),a(iptr3),a(iptr4),          &
       a(iptr5),a(iptr6),mxpts,nwhts,ibcvar,nbcvar,ndim)
!
!  restore tem space
!
  CALL resett
!
! go back to check if any more 2-d vector needs b.c. calculation
!
  GO TO 500
  600   CONTINUE
!
!------------------------------------------------------------------
! calculate b.c.'s for 3-d vectors
!------------------------------------------------------------------
!
  ndim=3
!
! initialize b.c. calculating status with value 1 for each vari.
!
  DO i=1,nxyz3d
    ibcsta(i)=1
  END DO
!
  700   CONTINUE
!
! search for variables that need b.c.'s and group them according to
! their staggering
!
  nbcvar=0
  DO ibc=1,nxyz3d
!
! look for 1st 3-d vector that needs b.c. calculation
!
    IF(ibdxyz(1,ibc) == 1.AND.ibdxyz(2,ibc) > 0.AND. ibcsta(ibc) == 1) THEN
      ibc2=ibdxyz(2,ibc)
      ibcvar(nbcvar+1)=ibc
      ibcvar(nbcvar+2)=ibc2
      nbcvar=nbcvar+2
      ibcsta(ibc )=0
      ibcsta(ibc2)=0
      EXIT
    END IF
  END DO
  92    CONTINUE
!
! if the b.c.'s for all 3-d vectors are calculated goto next step
!
  IF(nbcvar == 0)GO TO 800
!
! perform interpolations for 3-d vectors listed in ibcvar(nbcvar).
!
  mxpts=2*MAX( nx1+ny1,nx2+ny2 )
  iptr1=igetsp(4*mxpts)
  iptr2=igetsp(4*mxpts)
  iptr3=igetsp(4*mxpts)
  iptr4=igetsp(2*mxpts)
  iptr5=igetsp(nwhts*mxpts*2)
  iptr6=igetsp(4*mxpts)
!
  CALL itpbv2(mptr1,mptr2,a(iptr1),a(iptr2),a(iptr3),a(iptr4),          &
       a(iptr5),a(iptr6),mxpts,nwhts,ibcvar,nbcvar,ndim)
!
!  restore tem space
!
  CALL resett
!
! go back to check if any more 3-d vector needs b.c. calculation
!
  GO TO 700
  800   CONTINUE
  RETURN
  999   WRITE(6,'(''ERROR IN EXCHBC: THE LEVEL OF INPUT GRID'',         &
        & '' <1 OR >lfine, the grid was #'',i3)') mptr1
  RETURN
  998   WRITE(6,'(''ERROR IN EXCHBC: THE LEVEL OF INPUT GRID'',         &
        & '' <1 OR >lfine, the grid was #'',i3)') mptr2
  RETURN
  997   WRITE(6,'(''ERROR IN EXCHBC: THE INPUT GRIDS MUST ON SAME LEVEL'' &
        &  ,'' the levels were '',2I4)') level1, level2
  RETURN
END SUBROUTINE exchbc
!
!--------------------------------------------------------------------------
!

SUBROUTINE itpbs2(mptr1,mptr2,pts,irec,isrc,igsrc,whts,                 &
           pts1,mxpts,nwhts,ibcvar,nbcvar,ndim)
!
! exchange boundary values for scalar varibales in IBCVAR(NBCVAR)
! between two overlapping grids at the same level by interpolation.
!
  DIMENSION pts(2,mxpts),irec(2,mxpts,2),isrc(2,mxpts,2),               &
      igsrc(mxpts),whts(nwhts,mxpts,2),pts1(2,mxpts),ibcvar(nbcvar)
  INTEGER :: npts(2)
  LOGICAL :: change(2)

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'

  nw1=nint( SQRT(FLOAT(nwhts)) )
  nw2=nw1
!
! Calculate and store the coor's of p-points of the receive
! grid in array pts with their indecies stored in irec.
!
  IF(nbcvar == 0) RETURN
  CALL chkdim(ndim,'ITPBS2')

  ivar=ibcvar(1)

  DO igrid=1,2
    ig=igrid

    IF( igrid == 1) THEN
      mrec=mptr1
      msrc=mptr2
    ELSE
      mrec=mptr2
      msrc=mptr1
    END IF

    nx=node(5,mrec)
    ny=node(6,mrec)
    cosr=rnode(21,mrec)
    sinr=rnode(22,mrec)
    dxr=rnode(9 ,mrec)
    dyr=rnode(10,mrec)
    xor=rnode(1,mrec)
    yor=rnode(2,mrec)

    IF(ndim == 2) THEN
      xshift=stgxy(1,ivar)
      yshift=stgxy(2,ivar)
      nx1=nx-idmxy(1,ivar)
      ny1=ny-idmxy(2,ivar)
    ELSE
      xshift=stgxyz(1,ivar)
      yshift=stgxyz(2,ivar)
      nx1=nx-idmxyz(1,ivar)
      ny1=ny-idmxyz(2,ivar)
    END IF

    xc=xshift*dxr*cosr
    ys=yshift*dyr*sinr
    xs=xshift*dxr*sinr
    yc=yshift*dyr*cosr

    xorp=xor+xc-ys
    yorp=yor+xs+yc
!
! To calculate the absolute coordinates of the model boundary points
!
! Attention:
!
! Here we assume that in x direction the boundaries are at i=1 and nx1.
! and at j=1 and ny1 in y direction.
! If in a particular solver (model), the boundaries for certain
! variables are defined at points other the above, the use needs to
! change the following loop accordingly (i.e. change ist, jst, iend
! and jend). As a result, the part in routine UPDSCL that calculates
! the points for interior points updating also need change.
!
    ist=1
    jst=1
    iend=nx1
    jend=ny1
    DO j=jst,jend,jend-jst
      jj=(j-jst)/(jend-jst)*(iend-ist+1)
      DO i=ist,iend
        ii=jj+i-ist+1
        x=(i-1)*dxr
        y=(j-1)*dyr
        pts (1,ii) = xorp+x*cosr-y*sinr
        pts (2,ii) = yorp+x*sinr+y*cosr
        irec(1,ii,ig) = i
        irec(2,ii,ig) = j
      END DO
    END DO
    numpts=MAX(0, 2*(iend-ist+1))

    DO i=ist,iend,iend-ist
      ii=(i-ist)/(iend-ist)*(jend-jst-1)
      DO j=jst+1,jend-1
        jj=ii+j-jst+numpts
        x=(i-1)*dxr
        y=(j-1)*dyr
        pts (1,jj) = xorp+x*cosr-y*sinr
        pts (2,jj) = yorp+x*sinr+y*cosr
        irec(1,jj,ig) = i
        irec(2,jj,ig) = j
      END DO
    END DO
    numpts=numpts+MAX(0,2*(jend-jst-1))
!
! Find the source points for receive grid points pts
!
    ierfl=0
    iptcmp=1
    CALL getsrc(mrec,msrc,1,xshift,yshift,pts,numpts,                   &
                irec(1,1,ig),igsrc,ierfl,iptcmp)
!
    npts(ig)=numpts
    change(ig) = .true.
    IF(npts(ig) == 0) change(ig)=.false.

    IF(npts(ig) == 0) CYCLE
    DO ip=1,npts(ig)
      pts1(1,ip)=pts(1,ip)
      pts1(2,ip)=pts(2,ip)
    END DO
!
! Calculate weights of interpolation for points pts. The points are
! stored in pts1 so that pts remained unchanged after calwht.
!
    CALL calwht(mrec,msrc,pts1,isrc(1,1,ig),whts(1,1,ig),npts(ig),      &
                nwhts,ivar,ndim)

  END DO
!
! now we have receive points in irec, source points in isrc and interpolation
! weights wghts, for both interp from grid 2 to 1 and from grid 1 to 2.
! The actual interpolations are done in INTPSL for scalars listed in
! ibcvar. Interp bwteen the two grids are done at the same time
! so that packing unpacking only has to be done once.
!
  nx1=node(5 ,mptr1)
  ny1=node(6 ,mptr1)
  nz1=node(14,mptr1)
  nx2=node(5 ,mptr2)
  ny2=node(6 ,mptr2)
  nz2=node(14,mptr2)
  IF(nz1 /= nz2)GO TO 999
!
! Interpolate for the boundary values for variables in list ibcvar
!
  IF(npts(1) == 0.AND.npts(2) == 0) GO TO 500
  DO i=1,nbcvar
    ivar=ibcvar(i)
    IF(ndim == 2) THEN
      n3rd=1
      iptr1=igtnxy(mptr1,ivar,1)
      iptr2=igtnxy(mptr2,ivar,1)
    ELSE
      n3rd=nz1
      iptr1=igtxyz(mptr1,ivar,1)
      iptr2=igtxyz(mptr2,ivar,1)
    END IF

    IF(change(1)) CALL intrps(a(iptr1),a(iptr2),nx1,ny1,nx2,ny2,n3rd,   &
         irec(1,1,1),isrc(1,1,1),whts(1,1,1),npts(1),nw1,nw2)

    IF(change(2)) CALL intrps(a(iptr2),a(iptr1),nx2,ny2,nx1,ny1,n3rd,   &
         irec(1,1,2),isrc(1,1,2),whts(1,1,2),npts(2),nw1,nw2)
!
! if change is false, no need to pack the data back.
!
    IF(ndim == 2) THEN
      CALL retnxy(mptr1,ivar,1,iptr1,change(1))
      CALL retnxy(mptr2,ivar,1,iptr2,change(2))
    ELSE
      CALL retxyz(mptr1,ivar,1,iptr1,change(1))
      CALL retxyz(mptr2,ivar,1,iptr2,change(2))
    END IF
  END DO

  500   CONTINUE
  RETURN
  999   WRITE(6,'(''ERROR IN INTPS2: NZ FOR TWO GRIDS NOT EQUAL,'')')
  WRITE(6,'(''THEY WERE '',2I5,'', JOB ABORTED.'')')nz1,nz2
  STOP
END SUBROUTINE itpbs2
!
!--------------------------------------------------------------------------
!

SUBROUTINE itpbv2(mptr1,mptr2,pts,irec,isrc,igsrc,whts,                 &
           pts1,mxpts,nwhts,ibcvar,nbcvar,ndim)
!
! exchange boundary values for vectors in IBCVAR(NBCVAR)
! between two overlapping grids at the same level by interpolation.
!
! input list: mptr1, mptr2,mxpts,nwhts,ibcvar,nbcvar,ndim
! working arrays: pts,irec,isrc,igsrc,whts
! index implication for e.g. irec:
!
!   irec(2,     mxpts,      2     )
!     |        |         |
!  x or y   points   components
!
! similarly for isrc, igsrc,whts.
!
  DIMENSION pts(2,mxpts,2),irec(2,mxpts,2),isrc(2,mxpts,2),             &
            igsrc(mxpts,2),whts(nwhts,mxpts,2),pts1(2,mxpts,2),         &
            ibcvar(nbcvar),numpts(2),npts(2)

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  LOGICAL :: samstg, smstg2, change(2)

  nw1=nint( SQRT(FLOAT(nwhts)) )
  nw2=nw1
!
  IF(nbcvar == 0) RETURN
  CALL chkdim(ndim,'ITPBV2')
!
! To loop through all vectors, which has two components stored consecutively
! in nbcvar. The 1st compt. has to be the one in x direction.
! Here, interp weights are calculated for each vector.
!
  DO ivtr=1,nbcvar,2

    ivar1=ibcvar(ivtr)
    ivar2=ibcvar(ivtr+1)

    samstg=smstg2(ivar1,ivar2,ndim)
!
! find storage pointers for the receive grid variables (2 compnts of vector).
! uppack them if neccesary.
!
    IF(ndim == 2) THEN
      xshft = MAX( stgxy(1,ivar1),stgxy(1,ivar2) )
      yshft = MAX( stgxy(2,ivar1),stgxy(2,ivar2) )
      iptr11=igtnxy(mptr1,ivar1,1)
      iptr12=igtnxy(mptr1,ivar2,1)
      iptr21=igtnxy(mptr2,ivar1,1)
      iptr22=igtnxy(mptr2,ivar2,1)
    ELSE
      xshft = MAX( stgxyz(1,ivar1),stgxyz(1,ivar2) )
      yshft = MAX( stgxyz(2,ivar1),stgxyz(2,ivar2) )
      iptr11=igtxyz(mptr1,ivar1,1)
      iptr12=igtxyz(mptr1,ivar2,1)
      iptr21=igtxyz(mptr2,ivar1,1)
      iptr22=igtxyz(mptr2,ivar2,1)
    END IF

    DO igrid=1,2

      ig=igrid
      change(ig)=.true.
      IF(ig == 1)THEN
        mrec=mptr1
        msrc=mptr2
      ELSE
        mrec=mptr2
        msrc=mptr1
      END IF

      nxr=node(5,mrec)
      nyr=node(6,mrec)
      nzr=node(14,mrec)
      cosr=rnode(21,mrec)
      sinr=rnode(22,mrec)
      dxr=rnode(9 ,mrec)
      dyr=rnode(10,mrec)
      xor=rnode(1,mrec)
      yor=rnode(2,mrec)

      nxs=node(5,msrc)
      nys=node(6,msrc)
      nzs=node(14,msrc)
      coss=rnode(21,msrc)
      sins=rnode(22,msrc)

      cossr= coss*cosr+sins*sinr
      sinsr= coss*sinr-sins*cosr
      IF(nzr /= nzs) GO TO 999
      n3rd=nzs
      IF(ndim == 2) n3rd=1
!
! Find the points on rec. grid for both compnts. and the source grid for
! them.If the two compnts have the same staggering, they have the
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
            pts (1,ip,2)=pts (1,ip,1)
            pts (2,ip,2)=pts (2,ip,1)
            irec(1,ip,2)=irec(1,ip,1)
            irec(2,ip,2)=irec(2,ip,1)
            igsrc(ip,2)=igsrc(ip,1)
          END DO
          CYCLE
        END IF

        ivar=ibcvar(ivtr+icmpr-1)

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
! To calculate the absolute coordinates of the model boundary points
!
! Attention:
!
! Here we assume that in x direction the boundaries are at i=1 and nx1.
! and at j=1 and ny1 in y direction.
! If in a particular solver (model), the boundaries for certain
! variables are defined at points other the above, the use needs to
! change the following loop accordingly (i.e. change ist, jst, iend
! and jend). As a result, the part in routine UPDVTR that calculates
! the points for interior points updating also need change.
!
        ist=1
        jst=1
        iend=nx1
        jend=ny1
        DO j=jst,jend,jend-jst
          jj=(j-jst)/(jend-jst)*(iend-ist+1)
          DO i=ist,iend
            ii=jj+i-ist+1
            x=(i-1)*dxr
            y=(j-1)*dyr
            pts (1,ii,icmpr) = xorp+x*cosr-y*sinr
            pts (2,ii,icmpr) = yorp+x*sinr+y*cosr
            irec(1,ii,icmpr) = i
            irec(2,ii,icmpr) = j
          END DO
        END DO
        numpts(icmpr)=MAX(0,2*(iend-ist+1))

        DO i=ist,iend,iend-ist
          ii=(i-ist)/(iend-ist)*(jend-jst-1)
          DO j=jst+1,jend-1
            jj=ii+j-jst+numpts(icmpr)
            x=(i-1)*dxr
            y=(j-1)*dyr
            pts (1,jj,icmpr) = xorp+x*cosr-y*sinr
            pts (2,jj,icmpr) = yorp+x*sinr+y*cosr
            irec(1,jj,icmpr) = i
            irec(2,jj,icmpr) = j
          END DO
        END DO
        numpts(icmpr)=numpts(icmpr)+MAX(0,2*(jend-jst-1))
!
! Find the source grid for receive grid points pts
!
        ierfl=0
        iptcmp=1
        CALL getsrc(mrec,msrc,1,xshft,yshft,pts(1,1,icmpr),             &
            numpts(icmpr),irec(1,1,icmpr),igsrc(1,icmpr),ierfl,iptcmp)
        npts(icmpr)=numpts(icmpr)
        IF(npts(icmpr) == 0) change(ig)=.false.
!
      END DO

      DO icmps=1,2
        DO  icmpr=1,2

          IF(icmpr == 2.AND.samstg)THEN
! if the two components have the same staggering, simply make a copy for the
! 2nd component.
            npts(2)=npts(1)
            IF(npts(2) == 0) CYCLE
            DO ip=1,npts(2)
              pts (1,ip,2)=pts (1,ip,1)
              pts (2,ip,2)=pts (2,ip,1)
              isrc(1,ip,2)=isrc(1,ip,1)
              isrc(2,ip,2)=isrc(2,ip,1)
              irec(1,ip,2)=irec(1,ip,1)
              irec(2,ip,2)=irec(2,ip,1)
            END DO
            DO iw=1,nwhts
              DO ip=1,npts(2)
                whts(iw,ip,2)=whts(iw,ip,1)
              END DO
            END DO
            CYCLE
          END IF

          IF(npts(icmpr) == 0) CYCLE
          DO ip=1,npts(icmpr)
            pts1(1,ip,icmpr)=pts(1,ip,icmpr)
            pts1(2,ip,icmpr)=pts(2,ip,icmpr)
          END DO
!
! Calculate weights of interpolation for points pts. The points are
! stored in pts1 so that pts remained unchanged after calwht.
!
          ivsrc=ibcvar(ivtr+icmps-1)
          CALL calwht(mrec,msrc,pts1(1,1,icmpr),isrc(1,1,icmpr),        &
                      whts(1,1,icmpr),npts(icmpr),nwhts,ivsrc,ndim)
        END DO
        IF(npts(1) == 0.AND.npts(2) == 0) CYCLE
!
! find pointer for one of the components (icmps) of the source grid
! vector. If it is packed, do uppacking for this variable.
!
        ivsrc=ibcvar(ivtr+icmps-1)

        IF(icmps == 1)THEN
          proj1=cossr
          proj2=-sinsr
        ELSE
          proj1=sinsr
          proj2=cossr
        END IF

        IF(ig == 1)THEN
          irptr1=iptr11
          irptr2=iptr12
          isptr1=iptr21
          isptr2=iptr22
        ELSE
          irptr1=iptr21
          irptr2=iptr22
          isptr1=iptr11
          isptr2=iptr12
        END IF
        IF(icmps == 1) THEN
          isptr=isptr1
        ELSE
          isptr=isptr2
        END IF
        icmpr=1
        IF(npts(icmpr) /= 0) THEN
          CALL intrpv(a(irptr1),a(isptr),nxr,nyr,nxs,nys,n3rd,          &
              irec(1,1,icmpr),isrc(1,1,icmpr),whts(1,1,icmpr),          &
              npts(icmpr),nw1,nw2,icmps,proj1)
        END IF

        icmpr=2
        IF(npts(icmpr) /= 0) THEN
          CALL intrpv(a(irptr2),a(isptr),nxr,nyr,nxs,nys,n3rd,          &
              irec(1,1,icmpr),isrc(1,1,icmpr),whts(1,1,icmpr),          &
              npts(icmpr),nw1,nw2,icmps,proj2)
        END IF

      END DO
    END DO
    IF(ndim == 2) THEN
      CALL retnxy(mptr1,ivar1,1,iptr11,change(1))
      CALL retnxy(mptr1,ivar2,1,iptr12,change(1))
      CALL retnxy(mptr2,ivar1,1,iptr21,change(2))
      CALL retnxy(mptr2,ivar2,1,iptr22,change(2))
    ELSE
      CALL retxyz(mptr1,ivar1,1,iptr11,change(1))
      CALL retxyz(mptr1,ivar2,1,iptr12,change(1))
      CALL retxyz(mptr2,ivar1,1,iptr21,change(2))
      CALL retxyz(mptr2,ivar2,1,iptr22,change(2))
    END IF
  END DO
  RETURN
  999   WRITE(6,'(''ERROR IN INTPV2: NZ FOR TWO GRIDS NOT EQUAL,'')')
  WRITE(6,'(''THEY WERE '',2I5,'', JOB ABORTED.'')')nzr,nzs
  STOP
END SUBROUTINE itpbv2
