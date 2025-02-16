!
!--------------------------------------------------------------------------
!

SUBROUTINE inivar(mrec,msrc,lisvar,nvar,ndim)
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
  INTEGER :: lisvar(nvar)
  LOGICAL :: smstg2
!
!  interp all variables in lisvar(1-->nvar) for grid mrec from
!  grid msrc.  ndim is the dimension of the variables (either 2d or 3d)
!
  IF( intrpodr == 2 ) THEN
    nwhts=9
  ELSE
    nwhts=4
  END IF

  IF(nvar <= 0) THEN
    WRITE(6,'('' WARNING: NO VARIABLE WAS INITIALIZED IN INIVAR '',     &
    &   '' since nvar='',i4)') nvar
    RETURN
  END IF

  ldiff = ABS(node(4,mrec)-node(4,msrc))
  IF(ldiff > 1)THEN
    WRITE(6,'('' WARNING: RECEIVE AND SOURCE GRIDS DIFFER'',            &
    &   '' in level by'',i4,'', it should be =< 1 '')') ldiff
    RETURN
  END IF

  CALL chkdim(ndim,'INIVAR')

  nx=node(5,mrec)
  ny=node(6,mrec)
!
! check if the variable is a vector
!
  nchk=0
  IF(ndim == 2)THEN
    DO ivar=1,nvar
      IF (inixy(2,lisvar(ivar)) /= 0) nchk=nchk+1
    END DO
  ELSE
    DO ivar=1,nvar
      IF (inixyz(2,lisvar(ivar)) /= 0) nchk=nchk+1
    END DO
  END IF
  IF(nchk /= 0) THEN
    IF(nchk == nvar)GO TO 5000
    WRITE(6,'('' ERROR: THERE WAS A MIXTURE OF SCALARS AND'')')
    WRITE(6,'('' VECTORS PASSED INTO INIVAR IN ARRAY'',                 &
    &   '' lisvar, the list were:'')')
    WRITE(6,'(10x,2i10)') (i,lisvar(i),i=1,nvar)
    WRITE(6,'('' JOB STOPPED IN INIVAR!'')')
    STOP
  END IF
!
! if they are scalars, check if they have the same staggering
! if not, abort the job.
!
  nchk=0
  IF(ndim == 2)THEN
    DO i=2,nvar
      ivar=lisvar(i)
      IF(.NOT.smstg2(ivar,lisvar(1),2)) nchk=nchk+1
    END DO
  ELSE
    DO i=2,nvar
      ivar=lisvar(i)
      IF(.NOT.smstg2(ivar,lisvar(1),3)) nchk=nchk+1
    END DO
  END IF
  IF(nchk /= 0)THEN
    WRITE(6,'('' ERROR: LIST OF SCALARS DO NOT HAVE SAME '',            &
    &     ''staggering, job stopped in inivar.'')')
    STOP
  END IF
!
! perform interpolations for scalars listed in lisvar(nvar).
!
    mxpts=nx*ny
    iptr1=igetsp(2*mxpts)
    iptr2=igetsp(2*mxpts)
    iptr3=igetsp(2*mxpts)
    iptr4=igetsp(  mxpts)
    iptr5=igetsp(nwhts*mxpts)
    iptr6=igetsp(2*mxpts)
!
    CALL itpnsl(mrec,msrc,lisvar,nvar,ndim,                             &
         a(iptr1),a(iptr2),a(iptr3),                                    &
         a(iptr4),a(iptr5),a(iptr6),mxpts,nwhts)
!
    CALL resett
    GO TO 999
!
    5000  CONTINUE
!
! if in lisvar are vectors, make sure the matching components
! are stored consecutively.
!
    DO i=1,nvar-1,2
      IF(ndim == 2)THEN
        IF(lisvar(i+1) /= ABS(inixy (2,lisvar(i)))) GO TO 991
      ELSE
        IF(lisvar(i+1) /= ABS(inixyz(2,lisvar(i)))) GO TO 991
      END IF
    END DO
!
! Interpolate vectors from source grid MSRC to receive grid MREC.
!
    mxpts=nx*ny
    iptr1=igetsp(4*mxpts)
    iptr2=igetsp(4*mxpts)
    iptr3=igetsp(4*mxpts)
    iptr4=igetsp(2*mxpts)
    iptr5=igetsp(2*nwhts*mxpts)
    iptr6=igetsp(4*mxpts)
!
    DO i=1,nvar-1,2
      ivar1=lisvar(i)
      ivar2=lisvar(i+1)
!
! perform interpolations for a vector (ivar1, ivar2).
! Note that ivar1 must the component paralell to the x axis.
!
      CALL itpnvt(mrec,msrc,ivar1,ivar2,ndim,                           &
           a(iptr1),a(iptr2),a(iptr3),a(iptr4),a(iptr5),                &
           a(iptr6),mxpts,nwhts)
    END DO
!
    CALL resett
    999   RETURN
!
    991   WRITE(6,'('' ERROR: IMPROPER LIST OF VECTOR COMPONENTS '',    &
    &       ''given in array lisvar, they were:'')')
    WRITE(6,'(5x,3i10)')(i,lisvar(i),lisvar(i+1),i=1,nvar-1,2)
    WRITE(6,'('' JOB STOPPED IN INIVAR'')')
    STOP
END SUBROUTINE inivar
!
!--------------------------------------------------------------------------
!

SUBROUTINE itpnsl(mrec,msrc,lisvar,nvar,ndim,                           &
           pts,irec,isrc,igsrc,whts,pts1,mxpts,nwhts)
!
  DIMENSION pts(2,mxpts),irec(2,mxpts),isrc(2,mxpts),                   &
            igsrc(mxpts),whts(nwhts,mxpts),pts1(2,mxpts),               &
            lisvar(nvar)

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
  nw1=nint( SQRT(FLOAT(nwhts)) )
  nw2=nw1
!
  IF(nvar == 0) RETURN
  CALL chkdim(ndim,'INTNSL')

  nxr=node(5,mrec)
  nyr=node(6,mrec)
  nzr=node(14,mrec)

  nxs=node(5,msrc)
  nys=node(6,msrc)
  nzs=node(14,msrc)

  cosr=rnode(21,mrec)
  sinr=rnode(22,mrec)
  dxr=rnode(9 ,mrec)
  dyr=rnode(10,mrec)
  xor=rnode(1,mrec)
  yor=rnode(2,mrec)

  ivar=lisvar(1)
  IF(ndim == 2) THEN
    xshift=stgxy (1,ivar)
    yshift=stgxy (2,ivar)
    nx1=nxr-idmxy (1,ivar)
    ny1=nyr-idmxy (2,ivar)
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
  CALL getsrc(mrec,msrc, 1 ,xshift,yshift,pts,numpts,                   &
       irec,igsrc,ierfl,iptcmp)
!
  IF( numpts /= npts0 ) THEN
    WRITE(6,'(a,i5,a,i5,a,i5)')                                         &
        ' In INTNSL, ' ,npts0-numpts, ' points on grid ',mrec,          &
        ' got no source from source grid ',msrc
    WRITE(6,'(a,i1,a,i5)')                                              &
        ' It was for ',ndim,'-D scalar No. ',ivar
  END IF
!
! Calculate weights of interpolation for points pts.
!
  DO ip=1,numpts
    pts1(1,ip)=pts(1,ip)
    pts1(2,ip)=pts(2,ip)
  END DO
  CALL calwht(mrec,msrc,pts1,isrc,whts,numpts,nwhts,ivar,ndim)
!
! Interpolate for variables in list lisvar from source grid MSRC
! to receive grid MREC.
!
  DO i=1,nvar
    ivar=lisvar(i)
    IF(ndim == 2) THEN
      n3rd=1
      irptr=igtnxy(mrec,ivar,1)
      isptr=igtnxy(msrc,ivar,1)
    ELSE
      n3rd=nzs
      irptr=igtxyz(mrec,ivar,1)
      isptr=igtxyz(msrc,ivar,1)
    END IF

    CALL intrps(a(irptr),a(isptr),nxr,nyr,nxs,nys,n3rd,                 &
        irec,isrc,whts,numpts,nw1,nw2)

    IF(ndim == 2) THEN
      CALL retnxy(mrec,ivar,1,irptr,.true.)
      CALL retnxy(msrc,ivar,1,isptr,.false.)
    ELSE
      CALL retxyz(mrec,ivar,1,irptr,.true.)
      CALL retxyz(msrc,ivar,1,isptr,.false.)
    END IF
  END DO
!
  500   CONTINUE
  RETURN
END SUBROUTINE itpnsl
!
!--------------------------------------------------------------------------
!

SUBROUTINE itpnvt(mrec,msrc,ivar1a,ivar2a,ndim,                         &
           pts,irec,isrc,igsrc,whts,pts1,mxpts,nwhts)
!
  DIMENSION pts(2,mxpts,2),irec(2,mxpts,2),isrc(2,mxpts,2),             &
      igsrc(mxpts,2),whts(nwhts,mxpts,2),numpts(2),pts1(2,mxpts,2)
  DIMENSION lisvar(2)

  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
  LOGICAL :: samstg,smstg2

  nw1=nint( SQRT(FLOAT(nwhts)) )
  nw2=nw1
!
  CALL chkdim(ndim,'INTNVT')
!
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

  IF(ndim == 2)THEN
    i1= inixy (2,ivar1a)
    i2= inixy (2,ivar2a)
  ELSE
    i1= inixyz(2,ivar1a)
    i2= inixyz(2,ivar2a)
  END IF

  IF(i1*i2 >= 0) GO TO 991
  IF(i1 > 0.AND.i2 < 0) THEN
    ivar1 = ivar1a
    ivar2 = ivar2a
  END IF
  IF(i1 < 0.AND.i2 > 0)THEN
    WRITE(6,'('' WARNING: THE ORDER OF TWO COMPONENTS OF VACTOR''       &
    &   ,'' passed into intbvt incorrect!'')')
    WRITE(6,'('' THEIR ORDER IS CHANGED'')')
    ivar1 = ivar2a
    ivar2 = ivar1a
  END IF

  lisvar(1)=ivar1
  lisvar(2)=ivar2
!
! check to see if the two components of the vector have same staggering
!
  samstg=smstg2(ivar1,ivar2,ndim)
  IF(ndim == 2) THEN
    xshft = MAX( stgxy(1,ivar1),stgxy(1,ivar2) )
    yshft = MAX( stgxy(2,ivar1),stgxy(2,ivar2) )
  ELSE IF(ndim == 3)THEN
    xshft = MAX( stgxyz(1,ivar1),stgxyz(1,ivar2) )
    yshft = MAX( stgxyz(2,ivar1),stgxyz(2,ivar2) )
  END IF
!
! Find the points on rec. grid for both compnts. and their source grid.
! If the two compnts have the same staggering, the coordinate
! points and the source grid only need to be computed once.
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

    ivar=lisvar(icmpr)

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
        pts (1,ii,icmpr) = xorp+x*cosr-y*sinr
        pts (2,ii,icmpr) = yorp+x*sinr+y*cosr
        irec(1,ii,icmpr) = i
        irec(2,ii,icmpr) = j
      END DO
    END DO
    numpts(icmpr)=(jend-jst+1)*(iend-ist+1)
    npts0=numpts(icmpr)
!
! Find the source for receive grid points pts
!
    ierfl=0
    iptcmp=1
    CALL getsrc(mrec,msrc, 1 ,xshft,yshft,pts(1,1,icmpr),               &
        numpts(icmpr),irec(1,1,icmpr),igsrc(1,icmpr),ierfl,iptcmp)
!
    IF( (numpts(icmpr) /= npts0) .AND. verbose6 ) THEN
      WRITE(6,'('' WARNING: IN INTNVT,'',I4,'' POINTS ON GRID '',I4,    &
      &   '' got no source from source grid '',i4)')                    &
          npts0-numpts(icmpr), mrec, msrc
      WRITE(6,'('' IT WAS FOR '',I1,''-D VECTOR COMPNT. '',2I4)')       &
          ndim,ivar
    END IF
  END DO
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

  DO icmps=1,2
    DO  icmpr=1,2

      IF(icmpr == 2.AND.samstg)THEN
        DO ip=1,numpts(2)
          isrc(1,ip,2)=isrc (1,ip,1)
          isrc(2,ip,2)=isrc (2,ip,1)
        END DO
        DO iw=1,nwhts
          DO ip=1,numpts(2)
            whts(iw,ip,2)=whts(iw,ip,1)
          END DO
        END DO
        CYCLE
      END IF
!
! Calculate weights of interpolation for points pts.
!
      DO ip=1,numpts(icmpr)
        pts1(1,ip,icmpr)=pts(1,ip,icmpr)
        pts1(2,ip,icmpr)=pts(2,ip,icmpr)
      END DO
      ivsrc=lisvar(icmps)
      CALL calwht(mrec,msrc,pts1(1,1,icmpr),isrc(1,1,icmpr),            &
           whts(1,1,icmpr),numpts(icmpr),nwhts,ivsrc,ndim)
    END DO
!
! find pointer for one of the components (icmps) of the source grid
! vector. If it is packed, do uppacking for this variable.
!
    ivsrc=lisvar(icmps)
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
    IF(numpts(icmpr) /= 0) THEN
      CALL intrpv(a(irptr1),a(isptr),nxr,nyr,nxs,nys,n3rd,              &
          irec(1,1,icmpr),isrc(1,1,icmpr),whts(1,1,icmpr),              &
          numpts(icmpr),nw1,nw2,icmps,proj1)
    END IF

    icmpr=2
    IF(numpts(icmpr) /= 0) THEN
      CALL intrpv(a(irptr2),a(isptr),nxr,nyr,nxs,nys,n3rd,              &
          irec(1,1,icmpr),isrc(1,1,icmpr),whts(1,1,icmpr),              &
          numpts(icmpr),nw1,nw2,icmps,proj2)
    END IF

    IF(ndim == 2) THEN
      CALL retnxy(msrc,ivsrc,1,isptr,.false.)
    ELSE
      CALL retxyz(msrc,ivsrc,1,isptr,.false.)
    END IF
  END DO
!
  IF(ndim == 2) THEN
    CALL retnxy(mrec,ivar1,1,irptr1,.true.)
    CALL retnxy(mrec,ivar2,1,irptr2,.true.)
  ELSE
    CALL retxyz(mrec,ivar1,1,irptr1,.true.)
    CALL retxyz(mrec,ivar2,1,irptr2,.true.)
  END IF
  RETURN
  991   CONTINUE
  WRITE(6,'('' ERROR: TWO VARIABLES PASSED INTO INTBVT ARE'',           &
  &   '' two componoents of a vector!'')')
  WRITE(6,'('' JOB STOPPED IN INTBVT.'')')
  STOP
END SUBROUTINE itpnvt

SUBROUTINE chkdim(ndim,subnam)
  CHARACTER (LEN=*) :: subnam
  IF(ndim /= 2.AND.ndim /= 3)THEN
    WRITE(6,'('' ERROR IN '',A,                                         &
    &   '' : variable DIMENSION given wrong!'')')subnam
    WRITE(6,'(''IT SHOULD BE 2 OR 3, INPUT WAS'',I4)')ndim
    WRITE(6,'(''JOB STOPPED IN '',A)')subnam
    STOP
  END IF
  RETURN
END SUBROUTINE chkdim
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE INISFCTYPS               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initindex(mptr, mparent)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Fill soil and vegetation types into grid 'mptr' from its parent
!  grid mparent.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  04/23/1997
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!  mptr       Grid ID number
!  mparent    Grid mptr's parent grid ID number
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: mptr, mparent
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
  INTEGER :: nstyps
  PARAMETER ( nstyps = 4 )

  INTEGER :: isoiltyp1p, isoiltyp2p, isoiltyp3p, isoiltyp4p
  INTEGER :: istypfrct1p,istypfrct2p,istypfrct3p,istypfrct4p
  INTEGER :: isoiltyp1,  isoiltyp2,  isoiltyp3,  isoiltyp4
  INTEGER :: istypfrct1, istypfrct2, istypfrct3, istypfrct4
  INTEGER :: ivegtyp, ivegtypp

  INTEGER :: i,j, ij,ijs
  INTEGER :: nx,ny, nxp,nyp
  INTEGER :: is,js

  REAL :: hx,hy, hxp,hyp
  REAL :: xs0,ys0,xs,ys
  REAL :: xs0p,ys0p
!
!-----------------------------------------------------------------------
!
!  Integer function to get the pointers to the storage.
!
!-----------------------------------------------------------------------
!
  INTEGER :: igtnxy
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nx  = node(5,mptr)
  ny  = node(6,mptr)
  nxp = node(5,mparent)
  nyp = node(6,mparent)

  hx  = rnode(9 ,mptr)
  hy  = rnode(10,mptr)
  xs0 = rnode(1,mptr) + hx/2
  ys0 = rnode(2,mptr) + hy/2

  hxp  = rnode(9 ,mparent)
  hyp  = rnode(10,mparent)
  xs0p = rnode(1,mparent)
  ys0p = rnode(2,mparent)

  isoiltyp1  = igtnxy(mptr,id_soiltyp,  1)
  isoiltyp2  = igtnxy(mptr,id_soiltyp+1,1)
  isoiltyp3  = igtnxy(mptr,id_soiltyp+2,1)
  isoiltyp4  = igtnxy(mptr,id_soiltyp+3,1)

  istypfrct1 = igtnxy(mptr,id_stypfrct,  1)
  istypfrct2 = igtnxy(mptr,id_stypfrct+1,1)
  istypfrct3 = igtnxy(mptr,id_stypfrct+2,1)
  istypfrct4 = igtnxy(mptr,id_stypfrct+3,1)

  ivegtyp    = igtnxy(mptr,id_vegtyp,1)

  isoiltyp1p = igtnxy(mparent,id_soiltyp,  1)
  isoiltyp2p = igtnxy(mparent,id_soiltyp+1,1)
  isoiltyp3p = igtnxy(mparent,id_soiltyp+2,1)
  isoiltyp4p = igtnxy(mparent,id_soiltyp+3,1)

  istypfrct1p = igtnxy(mparent,id_stypfrct,  1)
  istypfrct2p = igtnxy(mparent,id_stypfrct+1,1)
  istypfrct3p = igtnxy(mparent,id_stypfrct+2,1)
  istypfrct4p = igtnxy(mparent,id_stypfrct+3,1)

  ivegtypp    = igtnxy(mparent,id_vegtyp,1)

  DO j=1,ny
    DO i=1,nx
      xs = xs0 + (i-1)*hx
      ys = ys0 + (j-1)*hy

      is = MAX( 1, MIN( nxp-1, INT((xs-xs0p)/hxp)+1 ) )
      js = MAX( 1, MIN( nyp-1, INT((ys-ys0p)/hyp)+1 ) )

      ij  = (j-1)*nx + i
      ijs = (js-1)*nxp + is

      a(isoiltyp1+ij-1) = a(isoiltyp1p+ijs-1)
      a(isoiltyp2+ij-1) = a(isoiltyp2p+ijs-1)
      a(isoiltyp3+ij-1) = a(isoiltyp3p+ijs-1)
      a(isoiltyp4+ij-1) = a(isoiltyp4p+ijs-1)

      a(istypfrct1+ij-1) = a(istypfrct1p+ijs-1)
      a(istypfrct2+ij-1) = a(istypfrct2p+ijs-1)
      a(istypfrct3+ij-1) = a(istypfrct3p+ijs-1)
      a(istypfrct4+ij-1) = a(istypfrct4p+ijs-1)

      a(ivegtyp+ij-1)  = a(ivegtypp+ijs-1)
    END DO
  END DO

  CALL retnxy(mptr,id_soiltyp, 4,isoiltyp1, .true.)
  CALL retnxy(mptr,id_stypfrct,4,istypfrct1,.true.)
  CALL retnxy(mptr,id_vegtyp,  1,ivegtyp,   .true.)

  RETURN

END SUBROUTINE initindex
