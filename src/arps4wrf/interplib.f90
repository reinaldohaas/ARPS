!#######################################################################

SUBROUTINE interp3d(sdata,mxs,mxe,mys,mye,nx,ny,nz,istagger,sx,sy,dx,dy,&
                    tx,ty,m1s,m1e,m2s,m2e,interp_method,                &
                    tdata,ibgn,iend,jbgn,jend,ksize,istatus)

  USE module_wrfgrid_constants

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: mxs,mxe,mys,mye,nx,ny,nz
  INTEGER, INTENT(IN)  :: istagger
  REAL,    INTENT(IN)  :: sdata(mxs:mxe,mys:mye,nz)
  REAL,    INTENT(IN)  :: sx(mxs:mxe), sy(mys:mye)
  REAL,    INTENT(IN)  :: dx, dy        ! source grid spacing
  INTEGER, INTENT(IN)  :: m1s,m1e,m2s,m2e
  REAL,    INTENT(IN)  :: tx(m1s:m1e,m2s:m2e), ty(m1s:m1e,m2s:m2e)
  INTEGER, INTENT(IN)  :: interp_method
  INTEGER, INTENT(IN)  :: ibgn,iend,jbgn,jend,ksize

  REAL,    INTENT(OUT) :: tdata(ibgn:iend,jbgn:jend,ksize)
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i, j, k

  INTEGER :: istart_x, iend_x, jstart_y, jend_y

  INTEGER :: ksrc, kfake_zone_bottom, isrc, jsrc, iorg, jorg

  INTEGER :: iorder
  INTEGER, ALLOCATABLE :: iloc(:,:), jloc(:,:)
                     ! WRF points(dest) index in ARPS grid (source)
  REAL,    ALLOCATABLE :: dxfld(:)
  REAL,    ALLOCATABLE :: dyfld(:)
  REAL,    ALLOCATABLE :: rdxfld(:)
  REAL,    ALLOCATABLE :: rdyfld(:)
  REAL,    ALLOCATABLE :: slopey(:,:,:)
  REAL,    ALLOCATABLE :: alphay(:,:,:)
  REAL,    ALLOCATABLE :: betay (:,:,:)

!  INTEGER  :: mnx, mny
!  REAL,    ALLOCATABLE :: tem1(:,:)

!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

  INTERFACE
    REAL FUNCTION interp_4pnt(mxs,mxe,mys,mye,istart_x,iend_x,jstart_y,jend_y, &
                              invar,inx,iny,rx,ry,dx,dy,istatus,maskarray,maskval)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: mxs, mxe, mys, mye
      INTEGER, INTENT(IN) :: istart_x, iend_x, jstart_y, jend_y
      REAL,    INTENT(IN) :: invar(mxs:mxe,mys:mye)
      REAL,    INTENT(IN) :: inx(mxs:mxe), iny(mys:mye)
      REAL,    INTENT(IN) :: rx, ry, dx, dy
      REAL,    INTENT(IN), OPTIONAL :: maskarray(mxs:mxe,mys:mye)
      REAL,    INTENT(IN), OPTIONAL :: maskval

      INTEGER, INTENT(OUT) :: istatus
    END FUNCTION interp_4pnt
  END INTERFACE

  REAL :: arpspntint2d

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  istart_x = mxs   ! source array with valid data
  iend_x   = mxe-1
  jstart_y = mys
  jend_y   = mye-1

  IF ( loc_x == 1 ) THEN
    istart_x = 1
  ELSE IF (loc_x == nproc_x ) THEN
    iend_x   = nx-1
  END IF

  IF ( loc_y == 1 ) THEN
    jstart_y = 1
  ELSE IF (loc_y == nproc_y ) THEN
    jend_y = ny-1
  END IF

  IF (istagger == 1) iend_x = iend_x + 1
  IF (istagger == 2) jend_y = jend_y + 1

  kfake_zone_bottom = 1
  IF (ksize == nz) kfake_zone_bottom = 0

!-----------------------------------------------------------------------

  SELECT CASE (interp_method)
  CASE (NEARNEIGHBOR)
    WRITE(6,'(1x,a,I2)') 'ERROR: interpolation method is still not implemented - ',interp_method
    CALL arpsstop('ERROR: not implemented interpolation method.',1)
  CASE (FOUR_PNT)
    DO k = 1,ksize
      ksrc = k + kfake_zone_bottom
      DO j = jbgn, jend
        DO i = ibgn, iend
          tdata(i,j,k) = interp_4pnt(mxs,mxe,mys,mye,istart_x,iend_x,jstart_y,jend_y, &
                         sdata(:,:,ksrc),sx,sy,tx(i,j),ty(i,j),dx,dy,istatus)
        END DO
      END DO
    END DO

  CASE (BILINEAR, QUADRATIC)
    ALLOCATE(iloc(ibgn:iend,jbgn:jend), STAT = istatus)
    ALLOCATE(jloc(ibgn:iend,jbgn:jend), STAT = istatus)
    ALLOCATE(dxfld (istart_x:iend_x),   STAT = istatus)
    ALLOCATE(rdxfld(istart_x:iend_x),   STAT = istatus)
    ALLOCATE(dyfld (jstart_y:jend_y),   STAT = istatus)
    ALLOCATE(rdyfld(jstart_y:jend_y),   STAT = istatus)

    ALLOCATE(slopey(istart_x:iend_x,jstart_y:jend_y,nz), STAT = istatus)
    ALLOCATE(alphay(istart_x:iend_x,jstart_y:jend_y,nz), STAT = istatus)
    ALLOCATE(betay (istart_x:iend_x,jstart_y:jend_y,nz), STAT = istatus)

    CALL setdestijloc(m1s,m1e,m2s,m2e,tx,ty,mxs,mxe,mys,mye,sx,sy,      &
                      ibgn,iend,jbgn,jend,iloc,jloc)

    CALL setsrcdxdy  (mxs,mxe,mys,mye,sx,sy,                            &
                      istart_x,iend_x,jstart_y,jend_y,                  &
                      dxfld,dyfld,rdxfld,rdyfld)

    CALL setsrcdrvy  (mxs,mxe,mys,mye,1,nz,                             &
                      istart_x,iend_x,jstart_y,jend_y,1,nz,             &
                      dyfld,rdyfld,sdata,                               &
                      slopey,alphay,betay)

!    IF (mp_opt > 0) THEn
!      mnx = iend_x-istart_x+1
!      mny = jend_y-jstart_y+1
!      ALLOCATE(tem1(mnx*mny,nz), STAT = istatus)
!
!write(0,*) 'here',istatus,myproc
!      CALL mpsendrecv2dew(slopey,mnx,mny,nz,0,0,1,tem1)
!write(0,*) 'here1',myproc
!      CALL mpsendrecv2dns(slopey,mnx,mny,nz,0,0,2,tem1)
!write(0,*) 'here2',myproc
!
!      CALL mpsendrecv2dew(alphay,mnx,mny,nz,0,0,1,tem1)
!      CALL mpsendrecv2dns(alphay,mnx,mny,nz,0,0,2,tem1)
!
!      CALL mpsendrecv2dew(betay, mnx,mny,nz,0,0,1,tem1)
!      CALL mpsendrecv2dns(betay, mnx,mny,nz,0,0,2,tem1)
!
!      DEALLOCATE(tem1)
!    END IF

    iorder = interp_method-BILINEAR+1

    DO k = 1,ksize
      ksrc = k + kfake_zone_bottom
      DO j = jbgn, jend
        DO i = ibgn, iend
          tdata(i,j,k) = arpspntint2d( mxs,mxe,mys,mye,                 &
                                istart_x,iend_x,jstart_y,jend_y,        &
                                iorder,sx,sy,tx(i,j),ty(i,j),           &
                                iloc(i,j),jloc(i,j),sdata(:,:,ksrc),    &
                                dxfld,dyfld,rdxfld,rdyfld,              &
                                slopey,alphay,betay)

        END DO
      END DO
    END DO

    DEALLOCATE(iloc,jloc,                   STAT = istatus)
    DEALLOCATE(dxfld,rdxfld,dyfld,rdyfld,   STAT = istatus)
    DEALLOCATE(slopey,alphay, betay,        STAT = istatus)

  CASE (SIXTEEN_PNT)
    WRITE(6,'(1x,a,I2)') 'ERROR: interpolation method is still not implemented - ',interp_method
    CALL arpsstop('ERROR: not implemented interpolation method.',1)

  CASE ( ASSIGN_DIRECT )
    DO k = 1,ksize
      ksrc = k + kfake_zone_bottom
      DO j = jbgn, jend
        jsrc = j-jbgn+2
        DO i = ibgn, iend
          isrc = i-ibgn+2
          tdata(i,j,k) = sdata(isrc,jsrc,ksrc)
        END DO
      END DO
    END DO

  CASE ( EXTRACT_SUBDOMAIN )

    CALL find_nearest_start(mxs,mxe,mys,mye,                            &
                     istart_x,iend_x,jstart_y,jend_y,                   &
                     sx,sy,tx(ibgn,jbgn),ty(ibgn,jbgn),dx,dy,           &
                     iorg,jorg,istatus)

    DO k = 1,ksize
      ksrc = k + kfake_zone_bottom
      DO j = jbgn, jend
        jsrc = j-jbgn+jorg
        DO i = ibgn, iend
          isrc = i-ibgn+iorg
          tdata(i,j,k) = sdata(isrc,jsrc,ksrc)
        END DO
      END DO
    END DO

  CASE DEFAULT
    WRITE(6,'(1x,a,I2)') 'ERROR: unknown interpolation method - ',interp_method
    CALL arpsstop('ERROR: unknown interpolation method.',1)
  END SELECT

  RETURN
END SUBROUTINE interp3d

!#######################################################################

SUBROUTINE interp3d_mask(sdata,mxs,mxe,mys,mye,nx,ny,nz,istagger,sx,sy,dx,dy,&
                    tx,ty,m1s,m1e,m2s,m2e,interp_method,                &
                    tdata,ibgn,iend,jbgn,jend,ksize,                    &
                    srcmask,destmask,mskval,istatus)

  USE module_wrfgrid_constants

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: mxs,mxe,mys,mye,nx,ny,nz
  INTEGER, INTENT(IN)  :: istagger
  REAL,    INTENT(IN)  :: sdata(mxs:mxe,mys:mye,nz)
  REAL,    INTENT(IN)  :: sx(mxs:mxe), sy(mys:mye)
  REAL,    INTENT(IN)  :: dx, dy        ! source grid spacing
  INTEGER, INTENT(IN)  :: m1s,m1e,m2s,m2e
  REAL,    INTENT(IN)  :: tx(m1s:m1e,m2s:m2e), ty(m1s:m1e,m2s:m2e)
  INTEGER, INTENT(IN)  :: interp_method
  INTEGER, INTENT(IN)  :: ibgn,iend,jbgn,jend,ksize

  REAL,    INTENT(IN)  :: srcmask(mxs:mxe,mys:mye)
  REAL,    INTENT(IN)  :: destmask(ibgn:iend,jbgn:jend)
           ! Source            Destination
           ! = 0, no value     = 0, water
           ! = 1, land         = 1, land
           ! = 2, water
  REAL,    INTENT(IN)  :: mskval        ! mask value, default value over water
  REAL,    INTENT(OUT) :: tdata(ibgn:iend,jbgn:jend,ksize)
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i, j, k

  INTEGER :: istart_x, iend_x, jstart_y, jend_y

  INTEGER :: ksrc, kfake_zone_bottom, isrc, jsrc, iorg, jorg

  INTEGER :: iorder
  INTEGER, ALLOCATABLE :: iloc(:,:), jloc(:,:)
                     ! WRF points(dest) index in ARPS grid (source)
  REAL,    ALLOCATABLE :: dxfld(:)
  REAL,    ALLOCATABLE :: dyfld(:)
  REAL,    ALLOCATABLE :: rdxfld(:)
  REAL,    ALLOCATABLE :: rdyfld(:)
  REAL,    ALLOCATABLE :: slopey(:,:,:)
  REAL,    ALLOCATABLE :: alphay(:,:,:)
  REAL,    ALLOCATABLE :: betay (:,:,:)

!  INTEGER  :: mnx, mny
!  REAL,    ALLOCATABLE :: tem1(:,:)

!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

  INTERFACE
    REAL FUNCTION interp_4pnt(mxs,mxe,mys,mye,istart_x,iend_x,jstart_y,jend_y, &
                              invar,inx,iny,rx,ry,dx,dy,istatus,maskarray,maskval)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: mxs, mxe, mys, mye
      INTEGER, INTENT(IN) :: istart_x, iend_x, jstart_y, jend_y
      REAL,    INTENT(IN) :: invar(mxs:mxe,mys:mye)
      REAL,    INTENT(IN) :: inx(mxs:mxe), iny(mys:mye)
      REAL,    INTENT(IN) :: rx, ry, dx, dy
      REAL,    INTENT(IN), OPTIONAL :: maskarray(mxs:mxe,mys:mye)
      REAL,    INTENT(IN), OPTIONAL :: maskval

      INTEGER, INTENT(OUT) :: istatus

    END FUNCTION interp_4pnt
  END INTERFACE

  REAL :: arpspntint2d

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  istart_x = mxs   ! source array with valid data
  iend_x   = mxe-1
  jstart_y = mys
  jend_y   = mye-1

  IF ( loc_x == 1 ) THEN
    istart_x = 1
  ELSE IF (loc_x == nproc_x ) THEN
    iend_x   = nx-1
  END IF

  IF ( loc_y == 1 ) THEN
    jstart_y = 1
  ELSE IF (loc_y == nproc_y ) THEN
    jend_y = ny-1
  END IF

  IF (istagger == 1) iend_x = iend_x + 1
  IF (istagger == 2) jend_y = jend_y + 1

  kfake_zone_bottom = 1
  IF (ksize == nz) kfake_zone_bottom = 0

!-----------------------------------------------------------------------

  SELECT CASE (interp_method)
  CASE (NEARNEIGHBOR)
    WRITE(6,'(1x,a,I2)') 'ERROR: interpolation method is still not implemented - ',interp_method
    CALL arpsstop('ERROR: not implemented interpolation method.',1)
  CASE (FOUR_PNT)
    DO k = 1,ksize
      ksrc = k + kfake_zone_bottom
      DO j = jbgn, jend
        DO i = ibgn, iend
          IF (destmask(i,j) > 0) THEN  ! land
            tdata(i,j,k) = interp_4pnt(mxs,mxe,mys,mye,istart_x,iend_x,jstart_y,jend_y, &
                           sdata(:,:,ksrc),sx,sy,tx(i,j),ty(i,j),dx,dy,istatus,srcmask,2.0)
          ELSE                         ! water
            tdata(i,j,k) = mskval
          END IF
        END DO
      END DO
    END DO

  CASE (BILINEAR, QUADRATIC)
    ALLOCATE(iloc(ibgn:iend,jbgn:jend), STAT = istatus)
    ALLOCATE(jloc(ibgn:iend,jbgn:jend), STAT = istatus)
    ALLOCATE(dxfld (istart_x:iend_x),   STAT = istatus)
    ALLOCATE(rdxfld(istart_x:iend_x),   STAT = istatus)
    ALLOCATE(dyfld (jstart_y:jend_y),   STAT = istatus)
    ALLOCATE(rdyfld(jstart_y:jend_y),   STAT = istatus)

    ALLOCATE(slopey(istart_x:iend_x,jstart_y:jend_y,nz), STAT = istatus)
    ALLOCATE(alphay(istart_x:iend_x,jstart_y:jend_y,nz), STAT = istatus)
    ALLOCATE(betay (istart_x:iend_x,jstart_y:jend_y,nz), STAT = istatus)

    CALL setdestijloc(m1s,m1e,m2s,m2e,tx,ty,mxs,mxe,mys,mye,sx,sy,      &
                      ibgn,iend,jbgn,jend,iloc,jloc)

    CALL setsrcdxdy  (mxs,mxe,mys,mye,sx,sy,                            &
                      istart_x,iend_x,jstart_y,jend_y,                  &
                      dxfld,dyfld,rdxfld,rdyfld)

    CALL setsrcdrvy  (mxs,mxe,mys,mye,1,nz,                             &
                      istart_x,iend_x,jstart_y,jend_y,1,nz,             &
                      dyfld,rdyfld,sdata,                               &
                      slopey,alphay,betay)

!    IF (mp_opt > 0) THEn
!      mnx = iend_x-istart_x+1
!      mny = jend_y-jstart_y+1
!      ALLOCATE(tem1(mnx*mny,nz), STAT = istatus)
!
!write(0,*) 'here',istatus,myproc
!      CALL mpsendrecv2dew(slopey,mnx,mny,nz,0,0,1,tem1)
!write(0,*) 'here1',myproc
!      CALL mpsendrecv2dns(slopey,mnx,mny,nz,0,0,2,tem1)
!write(0,*) 'here2',myproc
!
!      CALL mpsendrecv2dew(alphay,mnx,mny,nz,0,0,1,tem1)
!      CALL mpsendrecv2dns(alphay,mnx,mny,nz,0,0,2,tem1)
!
!      CALL mpsendrecv2dew(betay, mnx,mny,nz,0,0,1,tem1)
!      CALL mpsendrecv2dns(betay, mnx,mny,nz,0,0,2,tem1)
!
!      DEALLOCATE(tem1)
!    END IF

    iorder = interp_method-BILINEAR+1

    DO k = 1,ksize
      ksrc = k + kfake_zone_bottom
      DO j = jbgn, jend
        DO i = ibgn, iend
          IF (destmask(i,j) > 0) THEN  ! land
            tdata(i,j,k) = arpspntint2d( mxs,mxe,mys,mye,               &
                                istart_x,iend_x,jstart_y,jend_y,        &
                                iorder,sx,sy,tx(i,j),ty(i,j),           &
                                iloc(i,j),jloc(i,j),sdata(:,:,ksrc),    &
                                dxfld,dyfld,rdxfld,rdyfld,              &
                                slopey,alphay,betay)
          ELSE                         ! over ocean
            tdata(i,j,k) = mskval
          END IF

        END DO
      END DO
    END DO

    DEALLOCATE(iloc,jloc,                   STAT = istatus)
    DEALLOCATE(dxfld,rdxfld,dyfld,rdyfld,   STAT = istatus)
    DEALLOCATE(slopey,alphay, betay,        STAT = istatus)

  CASE (SIXTEEN_PNT)
    WRITE(6,'(1x,a,I2)') 'ERROR: interpolation method is still not implemented - ',interp_method
    CALL arpsstop('ERROR: not implemented interpolation method.',1)

  CASE ( ASSIGN_DIRECT )
    DO k = 1,ksize
      ksrc = k + kfake_zone_bottom
      DO j = jbgn, jend
        jsrc = j-jbgn+2
        DO i = ibgn, iend
          isrc = i-ibgn+2
          IF (destmask(i,j) > 0) THEN
            tdata(i,j,k) = sdata(isrc,jsrc,ksrc)
          ELSE
            tdata(i,j,k) = mskval
          END IF
        END DO
      END DO
    END DO

  CASE ( EXTRACT_SUBDOMAIN )

    CALL find_nearest_start(mxs,mxe,mys,mye,                            &
                     istart_x,iend_x,jstart_y,jend_y,                   &
                     sx,sy,tx(ibgn,jbgn),ty(ibgn,jbgn),dx,dy,           &
                     iorg,jorg,istatus)

    DO k = 1,ksize
      ksrc = k + kfake_zone_bottom
      DO j = jbgn, jend
        jsrc = j-jbgn+jorg
        DO i = ibgn, iend
          isrc = i-ibgn+iorg
          IF (destmask(i,j) > 0) THEN
            tdata(i,j,k) = sdata(isrc,jsrc,ksrc)
          ELSE
            tdata(i,j,k) = mskval
          END IF
        END DO
      END DO
    END DO

  CASE DEFAULT
    WRITE(6,'(1x,a,I2)') 'ERROR: unknown interpolation method - ',interp_method
    CALL arpsstop('ERROR: unknown interpolation method.',1)
  END SELECT

  RETURN
END SUBROUTINE interp3d_mask

!#######################################################################

SUBROUTINE interp2d(sdata,mxs,mxe,mys,mye,nx,ny,istagger,sx,sy,dx,dy,   &
                    tx,ty,m1s,m1e,m2s,m2e,interp_method,                &
                    tdata,ibgn,iend,jbgn,jend,istatus,landsea_mask)

  USE module_wrfgrid_constants

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: mxs,mxe,mys,mye,nx,ny
  INTEGER, INTENT(IN)  :: istagger
  REAL,    INTENT(IN)  :: sdata(mxs:mxe,mys:mye)
  REAL,    INTENT(IN)  :: sx(mxs:mxe), sy(mys:mye)
  REAL,    INTENT(IN)  :: dx, dy        ! source grid spacing
  INTEGER, INTENT(IN)  :: m1s,m1e,m2s,m2e
  REAL,    INTENT(IN)  :: tx(m1s:m1e,m2s:m2e), ty(m1s:m1e,m2s:m2e)
  INTEGER, INTENT(IN)  :: interp_method
  INTEGER, INTENT(IN)  :: ibgn,iend,jbgn,jend

  REAL,    INTENT(OUT) :: tdata(ibgn:iend,jbgn:jend)
  INTEGER, INTENT(OUT) :: istatus

  REAL,    INTENT(IN), OPTIONAL :: landsea_mask(ibgn:iend,jbgn:jend)

!-----------------------------------------------------------------------

  INTEGER :: i, j, k

  INTEGER :: istart_x, iend_x, jstart_y, jend_y

  INTEGER :: iorder, isrc, jsrc, iorg, jorg

  INTEGER, ALLOCATABLE :: iloc(:,:), jloc(:,:)
                     ! WRF points(dest) index in ARPS grid (source)
  REAL,    ALLOCATABLE :: dxfld(:)
  REAL,    ALLOCATABLE :: dyfld(:)
  REAL,    ALLOCATABLE :: rdxfld(:)
  REAL,    ALLOCATABLE :: rdyfld(:)
  REAL,    ALLOCATABLE :: slopey(:,:)
  REAL,    ALLOCATABLE :: alphay(:,:)
  REAL,    ALLOCATABLE :: betay (:,:)

!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

  INTERFACE
    REAL FUNCTION interp_4pnt(mxs,mxe,mys,mye,istart_x,iend_x,jstart_y,jend_y, &
                              invar,inx,iny,rx,ry,dx,dy,istatus,maskarray,maskval)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: mxs, mxe, mys, mye
      INTEGER, INTENT(IN) :: istart_x, iend_x, jstart_y, jend_y
      REAL,    INTENT(IN) :: invar(mxs:mxe,mys:mye)
      REAL,    INTENT(IN) :: inx(mxs:mxe), iny(mys:mye)
      REAL,    INTENT(IN) :: rx, ry, dx, dy
      REAL,    INTENT(IN), OPTIONAL :: maskarray(mxs:mxe,mys:mye)
      REAL,    INTENT(IN), OPTIONAL :: maskval

      INTEGER, INTENT(OUT) :: istatus

    END FUNCTION interp_4pnt
  END INTERFACE

  REAL :: interp_nearneighbor
  REAL :: arpspntint2d

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  istart_x = mxs   ! source array with valid data
  iend_x   = mxe-1
  jstart_y = mys
  jend_y   = mye-1

  IF ( loc_x == 1 ) THEN
    istart_x = 1
  ELSE IF (loc_x == nproc_x ) THEN
    iend_x   = nx-1
  END IF

  IF ( loc_y == 1 ) THEN
    jstart_y = 1
  ELSE IF (loc_y == nproc_y ) THEN
    jend_y = ny-1
  END IF

  IF (istagger == 1) iend_x = iend_x + 1
  IF (istagger == 2) jend_y = jend_y + 1

!-----------------------------------------------------------------------

  SELECT CASE (interp_method)
  CASE (NEARNEIGHBOR, NEARNEIGHBOR1, NEARNEIGHBOR2)
    DO j = jbgn, jend
      DO i = ibgn, iend
        tdata(i,j) = interp_nearneighbor(mxs,mxe,mys,mye,               &
                       istart_x,iend_x,jstart_y,jend_y,                 &
                       sdata(:,:),sx,sy,tx(i,j),ty(i,j),dx,dy,istatus)
      END DO
    END DO

    IF (interp_method == NEARNEIGHBOR1 ) THEN
      DO j = jbgn, jend
        DO i = ibgn, iend
          IF ( ABS(tdata(i,j)-13) < 1.0E-3 .OR.  &
               ABS(tdata(i,j)-12) < 1.0E-3) THEN  ! ARPS water or ice
            tdata(i,j) = 0.0            ! sea
          ELSE
            tdata(i,j) = 1.0            ! land
          END IF
        END DO
      END DO
    END IF

    IF (interp_method == NEARNEIGHBOR2 ) THEN
      DO j = jbgn, jend
        DO i = ibgn, iend
          IF ( ABS(tdata(i,j)-12) < 1.0E-3 .AND.     &  ! ARPS ice
               ABS(landsea_mask(i,j)) < 1.0E-3 ) THEN   ! it is indeed sea in the WRF grid
            tdata(i,j) = 1.0            ! ice flag
          ELSE
            tdata(i,j) = 0.0            ! not ice
          END IF
        END DO
      END DO
    END IF

  CASE (FOUR_PNT)
    DO j = jbgn, jend
      DO i = ibgn, iend
        tdata(i,j) = interp_4pnt(mxs,mxe,mys,mye,istart_x,iend_x,jstart_y,jend_y, &
                       sdata(:,:),sx,sy,tx(i,j),ty(i,j),dx,dy,istatus)
      END DO
    END DO

  CASE (BILINEAR, QUADRATIC)

    ALLOCATE(iloc(ibgn:iend,jbgn:jend), STAT = istatus)
    ALLOCATE(jloc(ibgn:iend,jbgn:jend), STAT = istatus)
    ALLOCATE(dxfld (istart_x:iend_x),   STAT = istatus)
    ALLOCATE(rdxfld(istart_x:iend_x),   STAT = istatus)
    ALLOCATE(dyfld (jstart_y:jend_y),   STAT = istatus)
    ALLOCATE(rdyfld(jstart_y:jend_y),   STAT = istatus)

    ALLOCATE(slopey(istart_x:iend_x,jstart_y:jend_y), STAT = istatus)
    ALLOCATE(alphay(istart_x:iend_x,jstart_y:jend_y), STAT = istatus)
    ALLOCATE(betay (istart_x:iend_x,jstart_y:jend_y), STAT = istatus)

    CALL setdestijloc(m1s,m1e,m2s,m2e,tx,ty,mxs,mxe,mys,mye,sx,sy,      &
                      ibgn,iend,jbgn,jend,iloc,jloc)

    CALL setsrcdxdy  (mxs,mxe,mys,mye,sx,sy,                            &
                      istart_x,iend_x,jstart_y,jend_y,                  &
                      dxfld,dyfld,rdxfld,rdyfld)

    CALL setsrcdrvy  (mxs,mxe,mys,mye,1,1,                              &
                      istart_x,iend_x,jstart_y,jend_y,1,1,              &
                      dyfld,rdyfld,sdata,                               &
                      slopey,alphay,betay)

    iorder = interp_method - BILINEAR + 1
    DO j = jbgn, jend
      DO i = ibgn, iend
        tdata(i,j) = arpspntint2d( mxs,mxe,mys,mye,                     &
                              istart_x,iend_x,jstart_y,jend_y,          &
                              iorder,sx,sy,tx(i,j),ty(i,j),             &
                              iloc(i,j),jloc(i,j),sdata(:,:),           &
                              dxfld,dyfld,rdxfld,rdyfld,                &
                              slopey,alphay,betay)

      END DO
    END DO

    DEALLOCATE(iloc,jloc,                   STAT = istatus)
    DEALLOCATE(dxfld,rdxfld,dyfld,rdyfld,   STAT = istatus)
    DEALLOCATE(slopey,alphay, betay,        STAT = istatus)

  CASE (ASSIGN_DIRECT)        ! Exact copy, no fake zone
      DO j = jbgn, jend
        jsrc = j-jbgn+2
        DO i = ibgn, iend
          isrc = i-ibgn+2
          tdata(i,j) = sdata(isrc,jsrc)
        END DO
      END DO

  CASE (EXTRACT_SUBDOMAIN)   ! The source array has already been extended. Extract from the extended array
                             ! It is similar to NEAREST_NEIGHBOR method
      CALL find_nearest_start(mxs,mxe,mys,mye,                          &
                       istart_x,iend_x,jstart_y,jend_y,                 &
                       sx,sy,tx(ibgn,jbgn),ty(ibgn,jbgn),dx,dy,         &
                       iorg,jorg,istatus)

      DO j = jbgn, jend
        jsrc = j - jbgn + jorg
        DO i = ibgn, iend
          isrc = i - ibgn + iorg
          tdata(i,j) = sdata(isrc,jsrc)
        END DO
      END DO

  CASE (SIXTEEN_PNT)
    WRITE(6,'(1x,a,I2)') 'ERROR: interpolation method is still not implemented - ',interp_method
    CALL arpsstop('ERROR: not implemented interpolation method.',1)
  CASE DEFAULT
    WRITE(6,'(1x,a,I2)') 'ERROR: unknown interpolation method - ',interp_method
    CALL arpsstop('ERROR: unknown interpolation method.',1)
  END SELECT

  RETURN
END SUBROUTINE interp2d

!#######################################################################

REAL FUNCTION interp_4pnt(mxs,mxe,mys,mye,istart_x,iend_x,jstart_y,jend_y, &
                          invar,inx,iny,rx,ry,dx,dy,istatus,maskarray,maskval)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mxs, mxe, mys, mye
  INTEGER, INTENT(IN) :: istart_x, iend_x, jstart_y, jend_y
  REAL,    INTENT(IN) :: invar(mxs:mxe,mys:mye)
  REAL,    INTENT(IN) :: inx(mxs:mxe), iny(mys:mye)
  REAL,    INTENT(IN) :: rx, ry, dx, dy
  REAL,    INTENT(IN), OPTIONAL :: maskarray(mxs:mxe,mys:mye)
  REAL,    INTENT(IN), OPTIONAL :: maskval

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ixf, ixc, jyf, jyc
  REAL    :: fxfy, fxcy, cxfy, cxcy
  REAL    :: dfx, dfy, dcx, dcy

!-----------------------------------------------------------------------
  REAL  :: interp_nearneighbor

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0
!
! Note that inx(1) is different in each processor and since
! rx, inx(1) as well dx are float number, the result may be
! minor different between MPI and non-mpi. That is the cause
! for not digit-by-digit identical.
!
  ixf = FLOOR  ( (rx-inx(1))/dx ) + 1
  ixc = CEILING( (rx-inx(1))/dx ) + 1
!  ixc = ixf + 1

  jyf = FLOOR  ( (ry-iny(1))/dy ) + 1
  jyc = CEILING( (ry-iny(1))/dy ) + 1

  ! First, make sure that the point is contained in the source array
  IF (ixf < istart_x .OR. ixc > iend_x  .OR.          &
      jyf < jstart_y .OR. jyc > jend_y) THEN

    istatus = -1
    IF ( ixc > iend_x .AND. rx < inx(iend_x)+dx*0.1) THEN
      ixc = ixf
      istatus = 0
    END IF
    IF ( jyc > jend_y .AND. ry < iny(jend_y)+dy*0.1) THEN
      jyc = jyf
      istatus = 0
    END IF

    IF (istatus < 0) THEN
      WRITE(6,'(a)') 'ERROR: destination point is outside of the source array in interp_4pnt.'
      WRITE(6,'(4(a,I4))') ' istart_x -- iend_x = ',istart_x,', ',iend_x, &
                          '; jstart_y -- jend_y = ',jstart_y,', ',jend_y
      WRITE(6,'(4(a,I4))') ' ixf, ixc           = ',ixf,', ',ixc,         &
                          '; jyf, jyc           = ',jyf,', ',jyc
      CALL arpsstop('ERROR: in interp_4pnt',1)
!    istatus = -1
!    RETURN
    END IF
  END IF

  dfx = (rx-inx(ixf))/dx
  dfy = (ry-iny(jyf))/dy
  dcx = (inx(ixc)-rx)/dx
  dcy = (iny(jyc)-ry)/dy

  fxfy = MAX(0.0, 1.0 - sqrt(dfx**2+dfy**2))
  fxcy = MAX(0.0, 1.0 - sqrt(dfx**2+dcy**2))
  cxfy = MAX(0.0, 1.0 - sqrt(dcx**2+dfy**2))
  cxcy = MAX(0.0, 1.0 - sqrt(dcx**2+dcy**2))

  IF ( present(maskarray) .AND. present(maskval) ) THEN
    IF ( ABS(maskarray(ixf,jyf) - maskval) < 0.01) fxfy = 0.0
    IF ( ABS(maskarray(ixf,jyc) - maskval) < 0.01) fxcy = 0.0
    IF ( ABS(maskarray(ixc,jyf) - maskval) < 0.01) cxfy = 0.0
    IF ( ABS(maskarray(ixc,jyc) - maskval) < 0.01) cxcy = 0.0
  END IF

  ! If all four points are missing, try the next interpolation method in the sequence
  IF ( (fxfy+fxcy+cxfy+cxcy) <= 0.0) THEN
    interp_4pnt = interp_nearneighbor(mxs,mxe,mys,mye,istart_x,iend_x,jstart_y,jend_y, &
                                      invar,inx,iny,rx,ry,dx,dy,istatus)
  ELSE
    interp_4pnt = ( fxfy*invar(ixf,jyf) + fxcy*invar(ixf,jyc) +         &
                  cxfy*invar(ixc,jyf) + cxcy*invar(ixc,jyc) ) /         &
                ( fxfy + fxcy + cxfy + cxcy )
  END IF

  RETURN
END FUNCTION interp_4pnt

!#######################################################################

REAL FUNCTION interp_nearneighbor(mxs,mxe,mys,mye,istart_x,iend_x,jstart_y,jend_y, &
                          invar,inx,iny,rx,ry,dx,dy,istatus)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mxs, mxe, mys, mye
  INTEGER, INTENT(IN) :: istart_x, iend_x, jstart_y, jend_y
  REAL,    INTENT(IN) :: invar(mxs:mxe,mys:mxe)
  REAL,    INTENT(IN) :: inx(mxs:mxe), iny(mys:mxe)
  REAL,    INTENT(IN) :: rx, ry, dx, dy
  REAL,    INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ixn, jyn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  ixn = NINT ( (rx-inx(1))/dx ) + 1
  jyn = NINT ( (ry-iny(1))/dy ) + 1

  ! First, make sure that the point is contained in the source array
  IF (ixn < istart_x .OR. ixn > iend_x .OR.          &
      jyn < jstart_y .OR. jyn > jend_y ) THEN
    WRITE(6,'(1x,a)') 'ERROR: destination point is outside of the'//    &
                      ' source array in interp_nearneighbor.'
    istatus = -1
    RETURN
  END IF

  interp_nearneighbor = invar(ixn,jyn)

  RETURN
END FUNCTION interp_nearneighbor
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION PNTINT2D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

REAL FUNCTION arpspntint2d(mxs,mxe,mys,mye,                             &
                           istart_x,iend_x,jstart_y,jend_y,             &
                           iorder,vx,vy,xpnt,ypnt,iloc,jloc,var,        &
                           dxfld,dyfld,rdxfld,rdyfld,                   &
                           slopey,alphay,betay)

!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: mxs,mxe,mys,mye
  INTEGER, INTENT(IN)  :: istart_x,iend_x,jstart_y, jend_y
  INTEGER, INTENT(IN)  :: iorder
  REAL,    INTENT(IN)  :: vx(mxs:mxe), vy(mys:mye)
  REAL,    INTENT(IN)  :: xpnt, ypnt
  INTEGER, INTENT(IN)  :: iloc,jloc
  REAL,    INTENT(IN)  :: var(mxs:mxe,mys:mye)
  REAL,    INTENT(IN)  :: dxfld(istart_x:iend_x)
  REAL,    INTENT(IN)  :: dyfld(jstart_y:jend_y)
  REAL,    INTENT(IN)  :: rdxfld(istart_x:iend_x)
  REAL,    INTENT(IN)  :: rdyfld(jstart_y:jend_y)
  REAL,    INTENT(IN)  :: slopey(istart_x:iend_x,jstart_y:jend_y)
  REAL,    INTENT(IN)  :: alphay(istart_x:iend_x,jstart_y:jend_y)
  REAL,    INTENT(IN)  :: betay(istart_x:iend_x,jstart_y:jend_y)

!-----------------------------------------------------------------------
  INTEGER :: ii, jj
  REAL    :: delx, dely
  REAL    :: alpha, beta, rtwodx
  REAL    :: varm1,var00,varp1
  REAL    :: varint

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (iorder == 1) THEN   ! Bilinear interpolation
    ii = MIN(MAX(iloc,istart_x),iend_x-1)
    jj = MIN(MAX(jloc,jstart_y),jend_y-1)
    delx = (xpnt - vx(ii))
    dely = (ypnt - vy(jj))
    varint = (1.0-delx*rdxfld(ii))*(var(ii,  jj)+slopey(ii,  jj)*dely)+ &
                 (delx*rdxfld(ii))*(var(ii+1,jj)+slopey(ii+1,jj)*dely)
  ELSE
    ii=MIN(MAX(iloc,(istart_x+1)),(iend_x-1))
    jj=MIN(MAX(jloc,(jstart_y+1)),(jend_y-1))
    delx = (xpnt - vx(ii))
    dely = (ypnt - vy(jj))
!
!-----------------------------------------------------------------------
!
!    Stencil is ii-1 to ii+1 and jj-1 to jj + 1
!
!    Interpolate in y.
!
!-----------------------------------------------------------------------
!
    varm1=(alphay(ii-1,jj)*dely+betay(ii-1,jj))*dely+var(ii-1,jj)
    var00=(alphay(ii  ,jj)*dely+betay(ii  ,jj))*dely+var(ii  ,jj)
    varp1=(alphay(ii+1,jj)*dely+betay(ii+1,jj))*dely+var(ii+1,jj)
!
!-----------------------------------------------------------------------
!
!    Interpolate intermediate results in x.
!
!-----------------------------------------------------------------------
!
    rtwodx= 1./(dxfld(ii-1)+dxfld(ii))
    alpha = ((varp1-var00)*rdxfld(ii)+(varm1-var00)*rdxfld(ii-1))*rtwodx
    beta  = (varp1-var00)*rdxfld(ii) - dxfld(ii)*alpha
    varint= (alpha*delx+beta)*delx+var00
  END IF

  arpspntint2d = varint

  RETURN
END FUNCTION arpspntint2d

!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SETIJLOC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE setdestijloc(m1s,m1e,m2s,m2e,x2d,y2d,                        &
                        mxs,mxe,mys,mye,x_ext,y_ext,                    &
                        ibgn,iend,jbgn,jend,iloc,jloc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Find i,j indices in verfication grid of each forecast point
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!  Modified from setijloc in file src/adas/intfield.f90.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: m1s,m1e,m2s,m2e
  INTEGER, INTENT(IN)  :: mxs,mxe,mys,mye
  REAL,    INTENT(IN)  :: x2d(m1s:m1e,m2s:m2e)
  REAL,    INTENT(IN)  :: y2d(m1s:m1e,m2s:m2e)
  REAL,    INTENT(IN)  :: x_ext(mxs:mxe)
  REAL,    INTENT(IN)  :: y_ext(mys:mye)
  INTEGER, INTENT(IN)  :: ibgn,iend,jbgn,jend
  INTEGER, INTENT(OUT) :: iloc(ibgn:iend,jbgn:jend)
  INTEGER, INTENT(OUT) :: jloc(ibgn:iend,jbgn:jend)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,i_ext,j_ext
  INTEGER :: imid,jmid
  REAL    :: xmid,ymid
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  imid=(mxe-mxs)/2 + mxs
  xmid=x_ext(imid)
  jmid=(mye-mys)/2 + mys
  ymid=y_ext(jmid)
!
  DO j = jbgn, jend
    DO i = ibgn, iend
      IF(x2d(i,j) < xmid) THEN
        DO i_ext=imid,mxs,-1
          IF(x_ext(i_ext) <= x2d(i,j)) EXIT
        END DO
        iloc(i,j)=i_ext
      ELSE
        DO i_ext=imid,mxe
          IF(x_ext(i_ext) >= x2d(i,j)) EXIT
        END DO
        iloc(i,j)=i_ext-1
      END IF
!
      IF(y2d(i,j) < ymid) THEN
        DO j_ext=jmid,mxs,-1
          IF(y_ext(j_ext) <= y2d(i,j)) EXIT
        END DO
        jloc(i,j)=j_ext
      ELSE
        DO j_ext=jmid,mxe
          IF(y_ext(j_ext) >= y2d(i,j)) EXIT
        END DO
        jloc(i,j)=j_ext-1
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE setdestijloc
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SETDXDY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setsrcdxdy(mxs,mxe,mys,mye,x1d,y1d,                          &
                      istart_x,iend_x,jstart_y,jend_y,                  &
                      dxfld,dyfld,rdxfld,rdyfld)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Calculate the local delta-x, delta-y and their inverses.
!    Precalculating these variables speeds up later calculations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: mxs,mxe,mys,mye
  REAL,    INTENT(IN)  :: x1d(mxs:mxe), y1d(mys:mye)
  INTEGER, INTENT(IN)  :: istart_x,iend_x,jstart_y,jend_y
  REAL,    INTENT(OUT) :: dxfld(istart_x:iend_x)
  REAL,    INTENT(OUT) :: dyfld(jstart_y:jend_y)
  REAL,    INTENT(OUT) :: rdxfld(istart_x:iend_x)
  REAL,    INTENT(OUT) :: rdyfld(jstart_y:jend_y)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO i=istart_x,iend_x-1
    dxfld(i)  = (x1d(i+1)-x1d(i))
    rdxfld(i) = 1./(x1d(i+1)-x1d(i))
  END DO
  dxfld (iend_x) = dxfld (iend_x-1)
  rdxfld(iend_x) = rdxfld(iend_x-1)
  DO j=jstart_y,jend_y-1
    dyfld(j)  = (y1d(j+1)-y1d(j))
    rdyfld(j) = 1./(y1d(j+1)-y1d(j))
  END DO
  dyfld (jend_y) = dyfld (jend_y-1)
  rdyfld(jend_y) = rdyfld(jend_y-1)
  RETURN
END SUBROUTINE setsrcdxdy
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SETDRVY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE setsrcdrvy(mxs,mxe,mys,mye,mzs,mze,                          &
                      ibeg,iend,jbeg,jend,kbeg,kend,                    &
                      dyfld,rdyfld,var,                                 &
                      slopey,alphay,betay)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Calculate the coefficients of interpolating polynomials
!    in the y-direction.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: mxs,mxe,mys,mye,mzs,mze
  INTEGER, INTENT(IN) :: ibeg,iend,jbeg,jend,kbeg,kend
  REAL,    INTENT(IN) :: dyfld(jbeg:jend)
  REAL,    INTENT(IN) :: rdyfld(jbeg:jend)
  REAL,    INTENT(IN) :: var(mxs:mxe,mys:mye,mzs:mze)
  REAL,    INTENT(OUT) :: slopey(ibeg:iend,jbeg:jend,kbeg:kend)
  REAL,    INTENT(OUT) :: alphay(ibeg:iend,jbeg:jend,kbeg:kend)
  REAL,    INTENT(OUT) :: betay(ibeg:iend,jbeg:jend,kbeg:kend)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL    :: rtwody
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbeg,kend
    DO j=jbeg+1,jend-1
      DO i=ibeg,iend
        slopey(i,j,k)=(var(i,j+1,k)-var(i,j,k))*rdyfld(j)
        rtwody=1./(dyfld(j-1)+dyfld(j))
        alphay(i,j,k)=((var(i,j+1,k)-var(i,j,k))*rdyfld(j) +            &
                       (var(i,j-1,k)-var(i,j,k))*rdyfld(j-1))*rtwody
        betay(i,j,k)=(var(i,j+1,k)-var(i,j,k))*rdyfld(j) -              &
                     dyfld(j)*alphay(i,j,k)
      END DO
    END DO
    slopey(:,jbeg,k) = 0.0
    slopey(:,jend,k) = 0.0
    alphay(:,jbeg,k) = 0.0
    alphay(:,jend,k) = 0.0
    betay (:,jbeg,k) = 0.0
    betay (:,jend,k) = 0.0
  END DO

  RETURN
END SUBROUTINE setsrcdrvy

!#######################################################################

SUBROUTINE find_nearest_start(mxs,mxe,mys,mye,istart_x,iend_x,jstart_y,jend_y, &
                              inx,iny,rx,ry,dx,dy,iorg,jorg,istatus)
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: mxs, mxe, mys, mye
  INTEGER, INTENT(IN)  :: istart_x, iend_x, jstart_y, jend_y
  REAL,    INTENT(IN)  :: inx(mxs:mxe), iny(mys:mxe)
  REAL,    INTENT(IN)  :: rx, ry, dx, dy
  INTEGER, INTENT(OUT) :: iorg, jorg
  REAL,    INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ixn, jyn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  iorg = NINT ( (rx-inx(1))/dx ) + 1
  jorg = NINT ( (ry-iny(1))/dy ) + 1

  ! First, make sure that the point is contained in the source array
  IF (iorg < istart_x .OR. iorg > iend_x .OR.          &
      jorg < jstart_y .OR. jorg > jend_y ) THEN
    WRITE(6,'(a)') 'ERROR: destination point is outside of the source array in find_nearest_start.'
    istatus = -1
    RETURN
  END IF

  RETURN
END SUBROUTINE find_nearest_start
