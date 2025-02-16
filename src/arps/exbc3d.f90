!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE EXTBDTINI                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE extbdtini(nx,ny,nz,                                          &
           u,v,w, ptprt,pprt,qv,qscalar,                                &
           ptbar,pbar,                                                  &
           u0exb, v0exb, w0exb,  pt0exb, pr0exb, qv0exb, qscalar0exb,   &
           udtexb,vdtexb,wdtexb,ptdtexb,prdtexb,qvdtexb,qscalardtexb,   &
           bcrlx,                                                       &
           uexbc,  vexbc, wexbc, ptexbc, prexbc, qvexbc, qscalarexbc,   &
           tem1,ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in predicted variables from the first available external data
!  sets to calculate the linear time tendencies of the external data
!  set.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  8/10/94
!
!  MODIFICATION HISTORY:
!
!  08/30/1995 (Yuhe Liu)
!  Changed the initial boundary arrays, for restart run, from the
!  model arrays to the time interplated external boundary arrays.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!  nx         Number of grid points in the x-direction (east/west)
!  ny         Number of grid points in the y-direction (north/south)
!  nz         Number of grid points in the vertical
!
!  nx,ny,nz   Number of grid points in x, y, and z dir.
!
!  u          u-velocity
!  v          v-velocity
!  w          w-velocity
!  ptprt      Potential temperature perturbation
!  pprt       Pressure perturbation
!  qv         Specific humidity
!  qc         Cloud water mixing ratio (kg/kg)
!  qr         Rain water mixing ratio (kg/kg)
!  qi         Cloud ice mixing ratio (kg/kg)
!  qs         Snow mixing ratio (kg/kg)
!  qh         Hail mixing ratio (kg/kg)
!
!  ptbar      Base state potential temperature
!  pbar       Base state pressure
!
!  OUTPUT:
!
!  uexbc      EXBC u  array
!  vexbc      EXBC v  array
!  wexbc      EXBC w  array
!  ptexbc     EXBC pt array
!  prexbc     EXBC pr array
!  qvexbc     EXBC qv array
!  qcexbc     EXBC qc array
!  qrexbc     EXBC qr array
!  qiexbc     EXBC qi array
!  qsexbc     EXBC qs array
!  qhexbc     EXBC qh array
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
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations and COMMON blocks.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz       ! Number of grid points in x, y, and z dir.

  REAL :: u(nx,ny,nz)       ! u-velocity
  REAL :: v(nx,ny,nz)       ! v-velocity
  REAL :: w(nx,ny,nz)       ! w-velocity
  REAL :: ptprt(nx,ny,nz)   ! Potential temperature perturbation
  REAL :: pprt(nx,ny,nz)    ! Pressure perturbation
  REAL :: qv(nx,ny,nz)      ! Specific humidity
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: ptbar(nx,ny,nz)   ! Base state potential temperature
  REAL :: pbar(nx,ny,nz)    ! Base state pressure

  REAL :: u0exb (nx,ny,nz)  ! External boundary u-velocity field
  REAL :: v0exb (nx,ny,nz)  ! External boundary v-velocity field
  REAL :: w0exb (nx,ny,nz)  ! External boundary w-velocity field
  REAL :: pt0exb(nx,ny,nz)  ! External boundary pt field
  REAL :: pr0exb(nx,ny,nz)  ! External boundary p field
  REAL :: qv0exb(nx,ny,nz)  ! External boundary qv field
  REAL :: qscalar0exb(nx,ny,nz,nscalar)

  REAL :: udtexb (nx,ny,nz) ! Time tendency of external boundary u
  REAL :: vdtexb (nx,ny,nz) ! Time tendency of external boundary v
  REAL :: wdtexb (nx,ny,nz) ! Time tendency of external boundary w
  REAL :: ptdtexb(nx,ny,nz) ! Time tendency of external boundary pt
  REAL :: prdtexb(nx,ny,nz) ! Time tendency of external boundary p
  REAL :: qvdtexb(nx,ny,nz) ! Time tendency of external boundary qv
  REAL :: qscalardtexb(nx,ny,nz,nscalar)

  REAL :: bcrlx (nx,ny)     ! EXBC relaxation coefficients

  REAL :: uexbc (nx,ny,nz)  ! EXBC u  array
  REAL :: vexbc (nx,ny,nz)  ! EXBC v  array
  REAL :: wexbc (nx,ny,nz)  ! EXBC w  array
  REAL :: ptexbc(nx,ny,nz)  ! EXBC pt array
  REAL :: prexbc(nx,ny,nz)  ! EXBC p  array
  REAL :: qvexbc(nx,ny,nz)  ! EXBC qv array
  REAL :: qscalarexbc(nx,ny,nz,nscalar)

  REAL :: tem1(nx,ny,nz)

  INTEGER :: abststrt

  CHARACTER ( LEN = 256 ) :: filename
  INTEGER                 :: lfname
  CHARACTER ( LEN = 15  ) :: ctime
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, n, nq
  INTEGER :: istat, ierr, ireturn
  INTEGER :: iyr, imon,idy, ihr, imin,isec
  INTEGER :: iebc,iwbc,jnbc,jsbc, idist
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Make sure the external boundary fields match and fit the model
!  domain.
!
!-----------------------------------------------------------------------
!
  n = 2*ngbrz
!
!-----------------------------------------------------------------------
!
!  Fill the boundary relaxation varibles with zero.
!
!-----------------------------------------------------------------------
!
  DO j = 1, ny
    DO i = 1, nx
      bcrlx(i,j) = 0.
    END DO
  END DO

  iebc=1
!     iwbc=nx-1
!     jnbc=ny-1
  iwbc=nx-1+(nproc_x-1)*(nx-3)
  jnbc=ny-1+(nproc_y-1)*(ny-3)
  jsbc=1
  IF (mp_opt > 0) THEN
  END IF

  DO j = 1, ny-1
    DO i = 1, nx-1
!         idist=min(i-iebc,iwbc-i,jnbc-j,j-jsbc)
      idist=MIN(i+(loc_x-1)*(nx-3)-iebc,                                &
                iwbc-i-(loc_x-1)*(nx-3),                                &
                jnbc-j-(loc_y-1)*(ny-3),                                &
                j+(loc_y-1)*(ny-3)-jsbc)
      IF(idist < ngbrz ) bcrlx(i,j) = 1./(1.+(REAL(idist)/brlxhw)**2)
    END DO
  END DO

  IF (myproc == 0) THEN
    WRITE (6,'(/a)') 'Boundary relaxation coefficients'
    WRITE (6,'(/a,a)')                                                  &
        ' j\\i    1       2       3       4       5       6       7',   &
        '       8       9'
    WRITE (6,'(i3,9f8.5)') (j,(bcrlx(i,j),i=1,9),j=1,ny)
  END IF

  CALL ctim2abss( year,month,day,hour,minute,second, abstinit )

  abststop  = abstinit + nint(tstop)
  abststrt  = abstinit + nint(tstart)

  CALL getbcfn( abststrt, exbcname, tinitebd, tintvebd,                 &
                filename, lfname, istat )

  IF ( istat /= 0 ) THEN
    ireturn = 1
    RETURN
  END IF

  IF ( initopt == 2 ) THEN
    RETURN
  END IF

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx
        u0exb(i,j,k) = u(i,j,k)
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny
      DO i = 1, nx-1
        v0exb(i,j,k) = v(i,j,k)
      END DO
    END DO
  END DO

  DO k = 1, nz
    DO j = 1, ny-1
      DO i = 1, nx-1
        w0exb(i,j,k) = w(i,j,k)
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        pt0exb(i,j,k) = ptprt(i,j,k)
        pr0exb(i,j,k) = pprt(i,j,k)
        qv0exb(i,j,k) = qv(i,j,k)
        DO nq = 1,nscalar
          qscalar0exb(i,j,k,nq) = qscalar(i,j,k,nq)
        END DO
      END DO
    END DO
  END DO

  abstfcst0 = abststrt

!  blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,readstride
    IF(myproc >= i.AND.myproc <= i+readstride-1)THEN

      IF(mp_opt > 0 .AND. readsplit(FINDX_B) > 0) THEN

        CALL readsplitexbc(nx,ny,nz,filename,lfname, ctime,             &
                    uexbc,vexbc,wexbc,ptexbc,prexbc,                    &
                    qvexbc,qscalarexbc,ierr)
      ELSE

        CALL readexbc(nx,ny,nz,filename,lfname, ctime,                  &
                    uexbc,vexbc,wexbc,ptexbc,prexbc,                    &
                    qvexbc,qscalarexbc,tem1,ierr)
      END IF

    END IF
    IF (mp_opt > 0) CALL mpbarrier
  END DO


  IF ( ierr == 1 ) THEN
    ireturn = 2
    RETURN
  ELSE IF ( ierr == 2 ) THEN
    ireturn = 3
    RETURN
  END IF

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        ptexbc(i,j,k) = ptexbc(i,j,k) - ptbar(i,j,k)
        prexbc(i,j,k) = prexbc(i,j,k) - pbar(i,j,k)
      END DO
    END DO
  END DO

  READ (ctime, '(i4,2i2,1x,3i2)')                                       &
       iyr,imon,idy,ihr,imin,isec

  CALL ctim2abss( iyr,imon,idy,ihr,imin,isec, abstfcst )

  DO k = 1, nz
    DO j = 1, ny
      DO i = 1, nx
        udtexb (i,j,k) = 0.0   ! Initialize the tendencies with zero.
        vdtexb (i,j,k) = 0.0
        wdtexb (i,j,k) = 0.0
        ptdtexb(i,j,k) = 0.0
        prdtexb(i,j,k) = 0.0
        qvdtexb(i,j,k) = 0.0
        DO nq = 1,nscalar
          qscalardtexb(i,j,k,nq) = 0.0
        END DO
      END DO
    END DO
  END DO

  CALL exbcdt(nx,ny,nz,                                                 &
              u0exb,v0exb,w0exb,pt0exb,pr0exb,qv0exb,qscalar0exb,       &
              udtexb,vdtexb,wdtexb,ptdtexb,prdtexb,qvdtexb,qscalardtexb,&
              uexbc,vexbc,wexbc,ptexbc,prexbc,qvexbc,qscalarexbc)

  RETURN
END SUBROUTINE extbdtini
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXTBDT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE extbdt(nx,ny,nz, ptbar,pbar, ireturn,                        &
           u0exb,v0exb,w0exb,pt0exb,pr0exb,qv0exb,qscalar0exb,          &
           udtexb,vdtexb,wdtexb,ptdtexb,prdtexb,qvdtexb,qscalardtexb,   &
           uexbc,vexbc,wexbc,ptexbc,prexbc,qvexbc,qscalarexbc,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in predicted variables from external boundary file to
!  calculate the linear time tendencies of the external data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!  8/10/94 (Yuhe Liu)
!  Split the initial call to calculate the time tendency into another
!  subroutine, EXTBDTINI.
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    uexbc    EXBC u  array
!    vexbc    EXBC v  array
!    wexbc    EXBC w  array
!    ptexbc   EXBC pt array
!    prexbc   EXBC pr array
!    qvexbc   EXBC qv array
!    qcexbc   EXBC qc array
!    qrexbc   EXBC qr array
!    qiexbc   EXBC qi array
!    qsexbc   EXBC qs array
!    qhexbc   EXBC qh array
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
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations and COMMON blocks.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz       ! Number of grid points in x, y, and z dir.

  REAL :: ptbar(nx,ny,nz)   ! Base state potential temperature
  REAL :: pbar(nx,ny,nz)    ! Base state pressure

  REAL :: u0exb (nx,ny,nz)  ! External boundary u-velocity field
  REAL :: v0exb (nx,ny,nz)  ! External boundary v-velocity field
  REAL :: w0exb (nx,ny,nz)  ! External boundary w-velocity field
  REAL :: pt0exb(nx,ny,nz)  ! External boundary pt field
  REAL :: pr0exb(nx,ny,nz)  ! External boundary p field
  REAL :: qv0exb(nx,ny,nz)  ! External boundary qv field
  REAL :: qscalar0exb(nx,ny,nz,nscalar)

  REAL :: udtexb (nx,ny,nz) ! Time tendency of external boundary u
  REAL :: vdtexb (nx,ny,nz) ! Time tendency of external boundary v
  REAL :: wdtexb (nx,ny,nz) ! Time tendency of external boundary w
  REAL :: ptdtexb(nx,ny,nz) ! Time tendency of external boundary pt
  REAL :: prdtexb(nx,ny,nz) ! Time tendency of external boundary p
  REAL :: qvdtexb(nx,ny,nz) ! Time tendency of external boundary qv
  REAL :: qscalardtexb(nx,ny,nz,nscalar)

  REAL :: uexbc (nx,ny,nz)  ! EXBC u  array
  REAL :: vexbc (nx,ny,nz)  ! EXBC v  array
  REAL :: wexbc (nx,ny,nz)  ! EXBC w  array
  REAL :: ptexbc(nx,ny,nz)  ! EXBC pt array
  REAL :: prexbc(nx,ny,nz)  ! EXBC p  array
  REAL :: qvexbc(nx,ny,nz)  ! EXBC qv array
  REAL :: qscalarexbc(nx,ny,nz,nscalar)

  REAL :: tem1(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields.
!
!-----------------------------------------------------------------------
!
  CHARACTER ( LEN = 256 ) :: filename
  INTEGER                 :: lfname
  CHARACTER ( LEN = 15  ) :: ctime
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: istat, ierr, ireturn
  INTEGER :: nq

  REAL    :: tinc
  INTEGER :: abstcur

  INTEGER :: iyr, imon,idy, ihr, imin,isec
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
  abstcur = abstinit + nint(curtim+dtbig)

  IF ( abstcur <= abstfcst ) THEN
    ireturn = 0
    RETURN
  END IF

  CALL getbcfn( abstcur, exbcname, tinitebd, tintvebd,                  &
                filename, lfname, istat )

  IF ( istat /= 0 ) THEN
    ireturn = 1
    RETURN
  END IF

  tinc = REAL( abstfcst - abstfcst0 )

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx
        u0exb(i,j,k)  = u0exb(i,j,k) + udtexb(i,j,k) * tinc
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny
      DO i = 1, nx-1
        v0exb(i,j,k)  = v0exb(i,j,k) + vdtexb(i,j,k) * tinc
      END DO
    END DO
  END DO

  DO k = 1, nz
    DO j = 1, ny-1
      DO i = 1, nx-1
        w0exb(i,j,k)  = w0exb(i,j,k) + wdtexb(i,j,k) * tinc
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        pt0exb(i,j,k) = pt0exb(i,j,k) + ptdtexb(i,j,k) * tinc
        pr0exb(i,j,k) = pr0exb(i,j,k) + prdtexb(i,j,k) * tinc
        qv0exb(i,j,k) = qv0exb(i,j,k) + qvdtexb(i,j,k) * tinc
        DO nq = 1,nscalar
          qscalar0exb(i,j,k,nq) = qscalar0exb(i,j,k,nq)                 &
                                + qscalardtexb(i,j,k,nq) * tinc
        END DO
      END DO
    END DO
  END DO

!  blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,readstride
    IF(myproc >= i.AND.myproc <= i+readstride-1)THEN

      IF(mp_opt > 0 .AND. readsplit(FINDX_B) > 0) THEN

        CALL readsplitexbc(nx,ny,nz,filename,lfname, ctime,             &
                    uexbc,vexbc,wexbc,ptexbc,prexbc,                    &
                    qvexbc,qscalarexbc,ierr)
      ELSE

        CALL readexbc(nx,ny,nz,filename,lfname, ctime,                  &
                    uexbc,vexbc,wexbc,ptexbc,prexbc,                    &
                    qvexbc,qscalarexbc,tem1,ierr)
      END IF

    END IF
    IF (mp_opt > 0) CALL mpbarrier
  END DO

  IF ( ierr == 1 ) THEN
    ireturn = 2
    RETURN
  ELSE IF ( ierr == 2 ) THEN
    ireturn = 3
    RETURN
  END IF

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        ptexbc(i,j,k) = ptexbc(i,j,k) - ptbar(i,j,k)
        prexbc(i,j,k) = prexbc(i,j,k) - pbar(i,j,k)
      END DO
    END DO
  END DO

  abstfcst0 = abstfcst

  READ (ctime, '(i4,2i2,1x,3i2)')                                       &
       iyr,imon,idy,ihr,imin,isec

  CALL ctim2abss( iyr,imon,idy,ihr,imin,isec, abstfcst )

  CALL exbcdt(nx,ny,nz,                                                 &
              u0exb,v0exb,w0exb,pt0exb,pr0exb,qv0exb,qscalar0exb,       &
              udtexb,vdtexb,wdtexb,ptdtexb,prdtexb,qvdtexb,qscalardtexb,&
              uexbc,vexbc,wexbc,ptexbc,prexbc,qvexbc,qscalarexbc)

  ireturn = 0

  RETURN
END SUBROUTINE extbdt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXBCDT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE exbcdt(nx,ny,nz,                                             &
           u0exb,v0exb,w0exb,pt0exb,pr0exb,qv0exb,qscalar0exb,          &
           udtexb,vdtexb,wdtexb,ptdtexb,prdtexb,qvdtexb,qscalardtexb,   &
           uexbc,vexbc,wexbc,ptexbc,prexbc,qvexbc,qscalarexbc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the linear time-dependent tendencies of external
!  boundary data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
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
  INCLUDE 'globcst.inc'
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations and COMMON blocks.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz       ! Number of grid points in x, y, and z dir.

  REAL :: u0exb (nx,ny,nz)  ! External boundary u-velocity field
  REAL :: v0exb (nx,ny,nz)  ! External boundary v-velocity field
  REAL :: w0exb (nx,ny,nz)  ! External boundary w-velocity field
  REAL :: pt0exb(nx,ny,nz)  ! External boundary pt field
  REAL :: pr0exb(nx,ny,nz)  ! External boundary p field
  REAL :: qv0exb(nx,ny,nz)  ! External boundary qv field
  REAL :: qscalar0exb(nx,ny,nz,nscalar)

  REAL :: udtexb (nx,ny,nz) ! Time tendency of external boundary u
  REAL :: vdtexb (nx,ny,nz) ! Time tendency of external boundary v
  REAL :: wdtexb (nx,ny,nz) ! Time tendency of external boundary w
  REAL :: ptdtexb(nx,ny,nz) ! Time tendency of external boundary pt
  REAL :: prdtexb(nx,ny,nz) ! Time tendency of external boundary p
  REAL :: qvdtexb(nx,ny,nz) ! Time tendency of external boundary qv
  REAL :: qscalardtexb(nx,ny,nz,nscalar)

  REAL :: uexbc (nx,ny,nz)  ! EXBC u  array
  REAL :: vexbc (nx,ny,nz)  ! EXBC v  array
  REAL :: wexbc (nx,ny,nz)  ! EXBC w  array
  REAL :: ptexbc(nx,ny,nz)  ! EXBC pt array
  REAL :: prexbc(nx,ny,nz)  ! EXBC p  array
  REAL :: qvexbc(nx,ny,nz)  ! EXBC qv array
  REAL :: qscalarexbc(nx,ny,nz,nscalar)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: nq

  REAL    :: tinc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL checkdims(nx,ny,nz, nxebc,nyebc,nzebc, 'EXBCDT')

  tinc = FLOAT( abstfcst - abstfcst0 )

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx
        udtexb(i,j,k) = (uexbc(i,j,k)-u0exb(i,j,k))/tinc
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny
      DO i = 1, nx-1
        vdtexb(i,j,k) = (vexbc(i,j,k)-v0exb(i,j,k))/tinc
      END DO
    END DO
  END DO

  DO k = 1, nz
    DO j = 1, ny-1
      DO i = 1, nx-1
        wdtexb(i,j,k) = (wexbc(i,j,k)-w0exb(i,j,k))/tinc
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        ptdtexb(i,j,k) = (ptexbc(i,j,k)-pt0exb(i,j,k))/tinc
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        prdtexb(i,j,k) = (prexbc(i,j,k)-pr0exb(i,j,k))/tinc
      END DO
    END DO
  END DO

  IF ( qvbcrd /= 0 ) THEN

    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          qvdtexb(i,j,k) = (qvexbc(i,j,k)-qv0exb(i,j,k))/tinc
        END DO
      END DO
    END DO

  END IF

  DO nq = 1,nscalar

    IF ( qscalarbcrd(nq) /= 0 ) THEN

      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            qscalardtexb(i,j,k,nq) = (qscalarexbc(i,j,k,nq)           &
                                        - qscalar0exb(i,j,k,nq))/tinc
          END DO
        END DO
      END DO

    END IF
  END DO

  RETURN
END SUBROUTINE exbcdt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXBCUV                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE exbcuv(nx,ny,nz, time, u,v, u0exb,v0exb,udtexb,vdtexb)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the external boundary conditions for u and v.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    time     The time at which the boundary conditions of u and v are set.
!
!    u        u-velocity
!    v        v-velocity
!
!  OUTPUT:
!
!    u        u-velocity
!    v        v-velocity
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz       ! Number of grid points in x, y, and z dir.

  REAL :: time              ! The time at which the boundary conditions
                            ! of u and v are set.

  REAL :: u(nx,ny,nz)       ! v-velocity
  REAL :: v(nx,ny,nz)       ! v-velocity

  REAL :: u0exb(nx,ny,nz)   ! External boundary u-velocity field
  REAL :: v0exb(nx,ny,nz)   ! External boundary v-velocity field

  REAL :: udtexb(nx,ny,nz)  ! Time tendency of external boundary u
  REAL :: vdtexb(nx,ny,nz)  ! Time tendency of external boundary v
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k

  REAL :: tinc
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL checkdims(nx,ny,nz, nxebc,nyebc,nzebc, 'EXBCUV')

  tinc = time - FLOAT(abstfcst0-abstinit)

  IF (loc_x == 1) THEN
    DO k = 1, nz-1
      DO j = 1, ny-1
        u( 1,j,k) = u0exb( 1,j,k) + udtexb( 1,j,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_x == nproc_x) THEN
    DO k = 1, nz-1
      DO j = 1, ny-1
        u(nx,j,k) = u0exb(nx,j,k) + udtexb(nx,j,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_y == 1) THEN
    DO k = 1, nz-1
      DO i = 1, nx
        u(i,   1,k) = u0exb(i,   1,k) + udtexb(i,   1,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_y == nproc_y) THEN
    DO k = 1, nz-1
      DO i = 1, nx
        u(i,ny-1,k) = u0exb(i,ny-1,k) + udtexb(i,ny-1,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_x == 1) THEN
    DO k = 1, nz-1
      DO j = 1, ny
        v(   1,j,k) = v0exb(   1,j,k) + vdtexb(   1,j,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_x == nproc_x) THEN
    DO k = 1, nz-1
      DO j = 1, ny
        v(nx-1,j,k) = v0exb(nx-1,j,k) + vdtexb(nx-1,j,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_y == 1) THEN
    DO k = 1, nz-1
      DO i = 1, nx-1
        v(i, 1,k) = v0exb(i, 1,k) + vdtexb(i, 1,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_y == nproc_y) THEN
    DO k = 1, nz-1
      DO i = 1, nx-1
        v(i,ny,k) = v0exb(i,ny,k) + vdtexb(i,ny,k) * tinc
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE exbcuv
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXBCW                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE exbcw( nx,ny,nz, time, w, w0exb,wdtexb )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the external boundary conditions for w.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    time     The time at which the boundary condition of w is set.
!
!    w        w-velocity
!
!  OUTPUT:
!
!    w        w-velocity
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz       ! Number of grid points in x, y, and z dir.

  REAL :: time              ! The time at which the boundary condition
                            ! of w is set.

  REAL :: w(nx,ny,nz)       ! w-velocity

  REAL :: w0exb(nx,ny,nz)   ! External boundary w-velocity field

  REAL :: wdtexb(nx,ny,nz)  ! Time tendency of external boundary w
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k

  REAL :: tinc
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL checkdims(nx,ny,nz, nxebc,nyebc,nzebc, 'EXBCW')

  tinc = time - FLOAT(abstfcst0-abstinit)

  IF (loc_x == 1) THEN
    DO k = 1, nz-1
      DO j = 1, ny-1
        w(   1,j,k) = w0exb(   1,j,k) + wdtexb(   1,j,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_x == nproc_x) THEN
    DO k = 1, nz-1
      DO j = 1, ny-1
        w(nx-1,j,k) = w0exb(nx-1,j,k) + wdtexb(nx-1,j,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_y == 1) THEN
    DO k = 1, nz-1
      DO i = 1, nx-1
        w(i,   1,k) = w0exb(i,   1,k) + wdtexb(i,   1,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_y == nproc_y) THEN
    DO k = 1, nz-1
      DO i = 1, nx-1
        w(i,ny-1,k) = w0exb(i,ny-1,k) + wdtexb(i,ny-1,k) * tinc
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE exbcw
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXBCPT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE exbcpt( nx,ny,nz, time, ptprt, pt0exb,ptdtexb )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the external boundary conditions for ptprt.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    time     The time at which the boundary condition of ptprt is set.
!
!    ptprt    Potential temperature perturbation
!
!  OUTPUT:
!
!    ptprt    Potential temperature perturbation
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz       ! Number of grid points in x, y, and z dir.

  REAL :: time              ! The time at which the boundary condition
                            ! of ptprt is set.
  REAL :: ptprt(nx,ny,nz)   ! Perturbation potential temperature.

  REAL :: pt0exb(nx,ny,nz)  ! External boundary pt field
  REAL :: ptdtexb(nx,ny,nz) ! Time tendency of external boundary pt
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: tinc
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL checkdims(nx,ny,nz, nxebc,nyebc,nzebc, 'EXBCPT')

  tinc = time - FLOAT(abstfcst0-abstinit)

  CALL exbcs( nx,ny,nz, ptprt, pt0exb, ptdtexb, tinc)

  RETURN
END SUBROUTINE exbcpt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXBCP                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE exbcp( nx,ny,nz, time , pprt, pr0exb,prdtexb )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the external boundary conditions for pprt.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    time     The time at which the boundary condition of pprt is set.
!
!    pprt     Pressure perturbation
!
!  OUTPUT:
!
!    pprt     Pressure perturbation
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz       ! Number of grid points in x, y, and z dir.

  REAL :: time              ! The time at which the boundary condition
                            ! of pprt is set.

  REAL :: pprt(nx,ny,nz)    ! Pressure perturbation

  REAL :: pr0exb(nx,ny,nz)  ! External boundary p field
  REAL :: prdtexb(nx,ny,nz) ! Time tendency of external boundary pt
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: tinc
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL checkdims(nx,ny,nz, nxebc,nyebc,nzebc, 'EXBCP')

  tinc = time - FLOAT(abstfcst0-abstinit)

  CALL exbcs( nx,ny,nz, pprt, pr0exb, prdtexb, tinc)

  RETURN
END SUBROUTINE exbcp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXBCQ                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE exbcq(nx,ny,nz, qflag, time, q, q0exb,qdtexb)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the external boundary conditions for q.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    qflag    Indicator of water/ice species
!
!    time     The time at which the boundary condition of q is set.
!
!    q        One of the water or ice species.
!
!  OUTPUT:
!
!    q        q updated at the lateral boundary
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz       ! Number of grid points in x, y, and z dir.
  INTEGER :: qflag          ! Indicator of water/ice species

  REAL :: time              ! The time at which the boundary condition
                            ! of q is set.

  REAL :: q(nx,ny,nz)       ! Water vapor mixing ratio

  REAL :: q0exb(nx,ny,nz)   ! External boundary qv field
  REAL :: qdtexb(nx,ny,nz)  ! Time tendency of external boundary qv
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: tinc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL checkdims(nx,ny,nz, nxebc,nyebc,nzebc, 'EXBCQ')

  tinc = time - FLOAT(abstfcst0-abstinit)

  IF (qflag == 0 .AND. qvbcrd == 1) THEN   ! qflag can be 0, but qscalarbcrd
                                           ! starts from 1
     CALL exbcs(nx,ny,nz,q,q0exb,qdtexb,tinc)
  ELSE IF (qscalarbcrd(qflag) > 0) THEN

      CALL exbcs(nx,ny,nz,q,q0exb,qdtexb,tinc)

  END IF

  RETURN
END SUBROUTINE exbcq
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXBCS                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE exbcs( nx,ny,nz, s, s0, sdt, tinc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the external boundary conditions for a scalar s.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue based on Yuhe Liu's EXBCP.
!  8/12/95
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    s        A scalar whose boundary values is to be set
!    s0       s at a past time
!    sdt      Time tendency of s
!    tinc     Time increment between currnet s and s0
!
!  OUTPUT:
!
!    s        Boundary values of s.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz ! Number of grid points in x, y, and z dir.

  REAL :: s  (nx,ny,nz)  ! A scalar
  REAL :: s0 (nx,ny,nz)  ! s at a past time
  REAL :: sdt(nx,ny,nz)  ! Time tendency of s
  REAL :: tinc                    ! Time increment between s and s0
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF (loc_x == 1) THEN
    DO k = 1, nz-1
      DO j = 1, ny-1
        s(   1,j,k) = s0(   1,j,k) + sdt(   1,j,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_x == nproc_x) THEN
    DO k = 1, nz-1
      DO j = 1, ny-1
        s(nx-1,j,k) = s0(nx-1,j,k) + sdt(nx-1,j,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_y == 1) THEN
    DO k = 1, nz-1
      DO i = 1, nx-1
        s(i,   1,k) = s0(i,   1,k) + sdt(i,   1,k) * tinc
      END DO
    END DO
  END IF

  IF (loc_y == nproc_y) THEN
    DO k = 1, nz-1
      DO i = 1, nx-1
        s(i,ny-1,k) = s0(i,ny-1,k) + sdt(i,ny-1,k) * tinc
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE exbcs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BRLXUVW                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE brlxuvw( nx,ny,nz, dtbig1,                                   &
           u,v,w,rhostr,                                                &
           uforce,vforce,wforce,                                        &
           u0exb,v0exb,w0exb, udtexb,vdtexb,wdtexb,bcrlx,               &
           tem1,tem2,tem3,tem4,tem5,tem6 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the boundary relaxation and computational mixing for u,
!  v, and w in the boundary zone.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Xue Ming & Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates (m/s).
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!  INPUT/OUTPUT :
!
!    uforce   forcing terms in u-momentum equation (kg/(m*s)**2).
!             uforce= uforce_others + uforce_boundary
!    vforce   forcing terms in v-momentum equation (kg/(m*s)**2).
!             vforce= vforce_others + vforce_boundary
!    wforce   forcing terms in w-momentum equation (kg/(m*s)**2).
!             wforce= wforce_others + wforce_boundary
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL :: dtbig1

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: uforce(nx,ny,nz)     ! forcing terms in u-momentum equ (kg/(m*s)**2)
                               ! uforce = uforce_others + uforce_boundary

  REAL :: vforce(nx,ny,nz)     ! forcing terms in v-momentum equ (kg/(m*s)**2)
                               ! vforce = vforce_others + vforce_boundary

  REAL :: wforce(nx,ny,nz)     ! forcing terms in w-momentum equ (kg/(m*s)**2)
                               ! wforce = wforce_others + wforce_boundary

  REAL :: u0exb (nx,ny,nz)     ! External boundary u-velocity field
  REAL :: v0exb (nx,ny,nz)     ! External boundary v-velocity field
  REAL :: w0exb (nx,ny,nz)     ! External boundary w-velocity field

  REAL :: udtexb (nx,ny,nz)    ! Time tendency of external boundary u
  REAL :: vdtexb (nx,ny,nz)    ! Time tendency of external boundary v
  REAL :: wdtexb (nx,ny,nz)    ! Time tendency of external boundary w

  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

  REAL :: tem1(nx,ny,nz)       ! Temporary array
  REAL :: tem2(nx,ny,nz)       ! Temporary array
  REAL :: tem3(nx,ny,nz)       ! Temporary array
  REAL :: tem4(nx,ny,nz)       ! Temporary array
  REAL :: tem5(nx,ny,nz)       ! Temporary array
  REAL :: tem6(nx,ny,nz)       ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: tinc, temb
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model execution
  INCLUDE 'grid.inc'          ! Grid parameters
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL checkdims(nx,ny,nz, nxebc,nyebc,nzebc, 'EXBCP')

  tinc = (curtim +dtbig-2.*dtbig1) - FLOAT(abstfcst0-abstinit)

  CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx
        temb = u0exb(i,j,k) + udtexb(i,j,k) * tinc
        tem4(i,j,k) = ( u(i,j,k) - temb ) * tem1(i,j,k)
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny
      DO i = 1, nx-1
        temb = v0exb(i,j,k) + vdtexb(i,j,k) * tinc
        tem5(i,j,k) = ( v(i,j,k) - temb ) * tem2(i,j,k)
      END DO
    END DO
  END DO

!   DO 30 k = 1, nz
!   DO 30 j = 1, ny-1
!   DO 30 i = 1, nx-1
!     temb = w0exb(i,j,k) + wdtexb(i,j,k) * tinc
!     tem6(i,j,k) = ( w(i,j,k) - temb ) * tem3(i,j,k)
!30    CONTINUE

  CALL difxx(tem4, nx,ny,nz, 2,nx-1,1,ny-1,1,nz-1, dx, tem1)
  CALL difyy(tem4, nx,ny,nz, 2,nx-1,2,ny-2,1,nz-1, dy, tem2)

  DO k = 1, nz-1
    DO j = 2, ny-2
      DO i = 2, nx-1
        uforce(i,j,k) = uforce(i,j,k)                                   &
                      - cbcdmp * .5 * (bcrlx(i-1,j)+bcrlx(i,j))         &
                               * tem4(i,j,k)                            &
                      + cbcmixh* ( 1.+2.*(bcrlx(i-1,j)+bcrlx(i,j)) )    &
                               * ( tem1(i,j,k) + tem2(i,j,k) )
      END DO
    END DO
  END DO

  CALL difxx(tem5, nx,ny,nz, 2,nx-2,2,ny-1,1,nz-1, dx, tem1)
  CALL difyy(tem5, nx,ny,nz, 1,nx-1,2,ny-1,1,nz-1, dy, tem2)

  DO k = 1, nz-1
    DO j = 2, ny-1
      DO i = 2, nx-2
        vforce(i,j,k) = vforce(i,j,k)                                   &
                      - cbcdmp * .5 * (bcrlx(i,j-1)+bcrlx(i,j))         &
                               * tem5(i,j,k)                            &
                      + cbcmixh* ( 1. + 2.*(bcrlx(i,j-1)+bcrlx(i,j)) )  &
                               * ( tem1(i,j,k) + tem2(i,j,k) )
      END DO
    END DO
  END DO

!  CALL difxx(tem6, nx,ny,nz, 2,nx-2,1,ny-1,1,nz-1, dx, tem1)
!  CALL difyy(tem6, nx,ny,nz, 1,nx-1,2,ny-2,1,nz-1, dy, tem2)
!
!  DO 300 k = 1, nz
!  DO 300 j = 2, ny-2
!  DO 300 i = 2, nx-2
!    wforce(i,j,k) = wforce(i,j,k)
!    :                - cbcdmp * bcrlx(i,j) * tem6(i,j,k)
!    :                + cbcmixh* ( 1. + 4.*bcrlx(i,j) )
!    :                         * ( tem1(i,j,k) + tem2(i,j,k) )
!300   CONTINUE
!
  RETURN
END SUBROUTINE brlxuvw
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BRLXPT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE brlxpt( nx,ny,nz, dtbig1, ptprt,rhostr, ptforce,             &
           pt0exb,ptdtexb, bcrlx,                                       &
           tem1, tem2, tem3,tem4 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the boundary relaxation and computational mixing for
!  potential temperature perturbation, ptprt, in the boundary zone.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Xue Ming & Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!  8/15/95 (M. Xue)
!  BRLXS is now called to calculated the boundary relaxation term
!
!-----------------------------------------------------------------------
!
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ptprt    Perturbation potential temperature (K)
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!  INPUT/OUTPUT :
!
!    ptforce  Source terms in potential temperature equation (kg/(m*s)**2).
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL :: dtbig1

  REAL :: ptprt  (nx,ny,nz)    ! Perturbation potential temperature (K)

  REAL :: rhostr (nx,ny,nz)    ! Base state density rhobar times j3.

  REAL :: ptforce(nx,ny,nz)    ! Forcing term in ptprt euation

  REAL :: pt0exb(nx,ny,nz)     ! External boundary pt field
  REAL :: ptdtexb(nx,ny,nz)    ! Time tendency of external boundary pt
  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

  REAL :: tem1(nx,ny,nz)       ! Temporary array
  REAL :: tem2(nx,ny,nz)       ! Temporary array
  REAL :: tem3(nx,ny,nz)       ! Temporary array
  REAL :: tem4(nx,ny,nz)       ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: tinc
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc' ! Global constants that control model execution
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL checkdims(nx,ny,nz, nxebc,nyebc,nzebc, 'BRLXPT')

  tinc = (curtim+dtbig-2.*dtbig1) - FLOAT(abstfcst0-abstinit)

  CALL brlxs(nx,ny,nz,ptprt,pt0exb,ptdtexb,bcrlx,tinc,rhostr,tem4,      &
             tem1,tem2,tem3)

  DO k = 1, nz-1
    DO j = 2, ny-2
      DO i = 2, nx-2
        ptforce(i,j,k) = ptforce(i,j,k) + tem4(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE brlxpt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BRLXP                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE brlxp( nx,ny,nz, dtbig1, pprt,rhostr, pforce,                &
           pr0exb,prdtexb, bcrlx,                                       &
           tem1, tem2, tem3,tem4)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the boundary relaxation and computational mixing for
!  pprt in the boundary zone.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Xue Ming & Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!  8/15/95 (M. Xue)
!  BRLXS is now called to calculated the boundary relaxation term
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates (m/s).
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!  INPUT/OUTPUT :
!
!    pforce   forcing terms in pressure equation (kg/(m*s)**2).
!             pforce = pforce_others + pforce_boundary
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL :: dtbig1

  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: pforce(nx,ny,nz)     ! forcing terms in pert. pressure (Pascal)
                               ! pforce = pforce_others + pforce_boundary

  REAL :: pr0exb(nx,ny,nz)     ! External boundary p field
  REAL :: prdtexb(nx,ny,nz)    ! Time tendency of external boundary p
  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

  REAL :: tem1(nx,ny,nz)       ! Temporary array
  REAL :: tem2(nx,ny,nz)       ! Temporary array
  REAL :: tem3(nx,ny,nz)       ! Temporary array
  REAL :: tem4(nx,ny,nz)       ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: tinc
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc' ! Global constants that control model execution
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL checkdims(nx,ny,nz, nxebc,nyebc,nzebc, 'BRLXP')

  tinc = (curtim+dtbig-2.*dtbig1) - FLOAT(abstfcst0-abstinit)

  CALL brlxs(nx,ny,nz,pprt,pr0exb,prdtexb,bcrlx,tinc,rhostr,tem4,       &
             tem1,tem2,tem3)

  DO k = 1, nz-1
    DO j = 2, ny-1
      DO i = 2, nx-1
        pforce(i,j,k) = pforce(i,j,k) + tem4(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE brlxp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BRLXQ                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE brlxq( nx,ny,nz, dtbig1, qflag, q,rhostr, qsrc,              &
                  qv0exb,q0exb,qvdtexb,qdtexb,                          &
                  bcrlx,tem1, tem2, tem3 ,tem4)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the boundary relaxation and computational mixing for
!  a water/ice variable q in the boundary zone.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Xue Ming & Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!  8/15/95 (M. Xue)
!  Generalized for any wter/ice variable q, by calling BRLXS.
!
!  9/15/97 (M. Xue and G. Bassett)
!  Added in subroutine BRLXQ checks like (qsbcrd.eq.1) in the IF
!  tests so that boundary relaxation or smoothing is NOT done when
!  variable is not found in the exbc file (when q*bcrd=0).
!  Previously it was ralaxing the variables to the initial state,
!  therefore initial storms in the boundary zone tend to stick.
!
!-----------------------------------------------------------------------
!
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    dtbig1   Big time step size (s)
!    qflag    Indicator for water/ice variable
!    q        One of the water or ice variable
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!  INPUT/OUTPUT :
!
!    qsrc     Source terms in water/ice water equation (kg/(m*s)**2).
!             qvsrc = qvsrc_others + qvsrc_boundary
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc' ! Global constants that control model execution
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL    :: dtbig1
  INTEGER :: qflag             ! Indicator for water/ice variable

  REAL :: q   (nx,ny,nz)       ! Water/ice variable

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: qsrc(nx,ny,nz)       ! Source term in q euation

  REAL :: qv0exb(nx,ny,nz)     ! External boundary qv field

  REAL :: q0exb(nx,ny,nz)

  REAL :: qvdtexb(nx,ny,nz)    ! Time tendency of external boundary qv

  REAL :: qdtexb(nx,ny,nz)

  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

  REAL :: tem1(nx,ny,nz)       ! Temporary array
  REAL :: tem2(nx,ny,nz)       ! Temporary array
  REAL :: tem3(nx,ny,nz)       ! Temporary array
  REAL :: tem4(nx,ny,nz)       ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL    :: tinc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  CALL checkdims(nx,ny,nz, nxebc,nyebc,nzebc, 'BRLXQ')

  tinc = (curtim+dtbig-2.*dtbig1) - FLOAT(abstfcst0-abstinit)

  DO k = 1, nz
    DO j = 1, ny
      DO i = 1, nx
        tem4(i,j,k) = 0.0
      END DO
    END DO
  END DO

  IF( qvbcrd == 1 .AND. qflag == 0 ) THEN

    CALL brlxs(nx,ny,nz,q,qv0exb,qvdtexb,bcrlx,tinc,rhostr,tem4,        &
               tem1,tem2,tem3)

  ELSE IF( qscalarbcrd(qflag) > 0 ) THEN

    CALL brlxs(nx,ny,nz,q,q0exb,qdtexb,bcrlx,tinc,rhostr,tem4,          &
               tem1,tem2,tem3)

  END IF

  DO k = 1, nz-1
    DO j = 2, ny-2
      DO i = 2, nx-2
        qsrc(i,j,k) = qsrc(i,j,k) + tem4(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE brlxq
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BRLXS                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE brlxs( nx,ny,nz, s, s0,sdt,bcrlx, tinc, rhostr, srlx,        &
           tem1, tem2, tem3 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the boundary relaxation and computational mixing term
!  for a scalar 's' in the boundary zone.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  8/15/95
!
!  MODIFICATION HISTORY:
!
!  8/24/95 (M. Xue)
!  Corrected a sign error in loop 60.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    s        A scalar whose boundary values is to be set
!    s0       s at a past time
!    sdt      Time tendency of s
!    tinc     Time increment between currnet s and s0
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!  OUTPUT:
!
!    srlx     Relaxation and spatial smoothing term in the boundary zone
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: s  (nx,ny,nz)        ! A scalar
  REAL :: s0 (nx,ny,nz)        ! s at a past time
  REAL :: sdt(nx,ny,nz)        ! Time tendency of s
  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

  REAL :: tinc                 ! Time increment between s and s0
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: srlx(nx,ny,nz)       ! Source term in qv euation

  REAL :: tem1(nx,ny,nz)       ! Temporary array
  REAL :: tem2(nx,ny,nz)       ! Temporary array
  REAL :: tem3(nx,ny,nz)       ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc' ! Global constants that control model execution
  INCLUDE 'grid.inc'          ! Grid parameters
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        tem1(i,j,k)=(s(i,j,k)-(s0(i,j,k)+sdt(i,j,k)*tinc))*rhostr(i,j,k)
      END DO
    END DO
  END DO

  CALL difxx(tem1, nx,ny,nz, 2,nx-2,1,ny-1,2,nz-1, dx, tem2)
  CALL difyy(tem1, nx,ny,nz, 1,nx-1,2,ny-2,2,nz-1, dy, tem3)

  DO k = 2, nz-2
    DO j = 2, ny-2
      DO i = 2, nx-2
        srlx(i,j,k) = - cbcdmp * bcrlx(i,j) * tem1(i,j,k)               &
                      + cbcmixh* ( 1. + 4.*bcrlx(i,j) )                 &
                               * ( tem2(i,j,k) + tem3(i,j,k) )
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE brlxs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE CHECKDIMS                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkdims(nx,ny,nz, nxebc,nyebc,nzebc, subname)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check if two sets of grid dimensions are the same.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  8/12/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!  nx         Number of grid points in the x-direction (east/west)
!  ny         Number of grid points in the y-direction (north/south)
!  nz         Number of grid points in the vertical
!  nxebc      nx in the external boundary condition data
!  nyebc      ny in the external boundary condition data
!  nzebc      nz in the external boundary condition data
!
!  OUTPUT:
!
!  None
!
!-----------------------------------------------------------------------
!
!  Variable Declarations and COMMON blocks.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  INTEGER :: nxebc,nyebc,nzebc
  CHARACTER (LEN=*) :: subname
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (nxebc /= nx .OR. nyebc /= ny .OR. nzebc /= nz) THEN
    WRITE (6,'(a/a/a/a,i5,a,i5,a,i5/a/a,i5,a,i5,a,i5)')                 &
        ' Array dimension(s) of the external boundary fields ',         &
        ' inconsistent with model definitions. ',                       &
        ' Dimensions for boundary fields were',                         &
        ' nxebc = ', nxebc, ', nyebc = ', nyebc,                        &
        ', nzebc = ', nzebc,                                            &
        ' and the model definitions were',                              &
        ' nx  = ', nx,  ', ny  = ', ny,  ', nz  = ', nz
    CALL arpsstop('arpsstop called from checkdims',1)
  END IF

  RETURN
END SUBROUTINE checkdims
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE SETEXBCCST               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setexbcptr(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the pointers to EXBC variables in the EXBC buffer array.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  04/17/1997
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz       ! The number grid points in 3 directions
  INTEGER :: nq

  INTEGER :: nxy, nxyz
  INTEGER :: nptr
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'exbc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxy  = nx*ny
  nxyz = nxy*nz
!
!-----------------------------------------------------------------------
!
!  Pointers to previous EXBC data arrays
!
!-----------------------------------------------------------------------
!
  nu0exb  = 1 + 00*nxyz
  nv0exb  = 1 + 01*nxyz
  nw0exb  = 1 + 02*nxyz
  npt0exb = 1 + 03*nxyz
  npr0exb = 1 + 04*nxyz
  nqv0exb = 1 + 05*nxyz

  nqscalar0exb(:) = 0
  nqscalar0exb(0) = nqv0exb

  DO nq = 1,nscalar
    nqscalar0exb(nq) = nqv0exb + nq*nxyz
  END DO

!  nqc0exb = 0
!  nqr0exb = 0
!  nqi0exb = 0
!  nqs0exb = 0
!  nqh0exb = 0
!  IF (P_QC > 0) nqc0exb = nqscalar0exb(P_QC)
!  IF (P_QR > 0) nqr0exb = nqscalar0exb(P_QR)
!  IF (P_QI > 0) nqi0exb = nqscalar0exb(P_QI)
!  IF (P_QS > 0) nqs0exb = nqscalar0exb(P_QS)
!  IF (P_QH > 0) nqh0exb = nqscalar0exb(P_QH)
!
!-----------------------------------------------------------------------
!
!  Pointers to previous EXBC time tendency arrays
!
!-----------------------------------------------------------------------
!

  nptr     = 06 + nscalar

  nudtexb  = 1 +     nptr*nxyz
  nvdtexb  = 1 + (nptr+1)*nxyz
  nwdtexb  = 1 + (nptr+2)*nxyz
  nptdtexb = 1 + (nptr+3)*nxyz
  nprdtexb = 1 + (nptr+4)*nxyz
  nqvdtexb = 1 + (nptr+5)*nxyz

  nqscalardtexb(:) = 0
  nqscalardtexb(0) = nqvdtexb

  DO nq = 1, nscalar
    nqscalardtexb(nq) = nqvdtexb + nq*nxyz
  END DO

!  nqcdtexb = 0
!  nqrdtexb = 0
!  nqidtexb = 0
!  nqsdtexb = 0
!  nqhdtexb = 0
!  IF (P_QC > 0) nqcdtexb = nqscalardtexb(P_QC)
!  IF (P_QR > 0) nqrdtexb = nqscalardtexb(P_QR)
!  IF (P_QI > 0) nqidtexb = nqscalardtexb(P_QI)
!  IF (P_QS > 0) nqsdtexb = nqscalardtexb(P_QS)
!  IF (P_QH > 0) nqhdtexb = nqscalardtexb(P_QH)

  RETURN
END SUBROUTINE setexbcptr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BRLXUVW_RBC                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE brlxuvw_rbc( nx,ny,nz, u,v,w,ubar,vbar,rhostr,               &
           uforce,vforce,wforce, bcrlx)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the boundary relaxation and computational mixing for u,
!  v, and w in the boundary zone.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Xue Ming & Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates (m/s).
!
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!  INPUT/OUTPUT :
!
!    uforce   forcing terms in u-momentum equation (kg/(m*s)**2).
!             uforce= uforce_others + uforce_boundary
!    vforce   forcing terms in v-momentum equation (kg/(m*s)**2).
!             vforce= vforce_others + vforce_boundary
!    wforce   forcing terms in w-momentum equation (kg/(m*s)**2).
!             wforce= wforce_others + wforce_boundary
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ubar  (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Total v-velocity (m/s)

  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.

  REAL :: uforce(nx,ny,nz)     ! forcing terms in u-momentum equ (kg/(m*s)**2)
                               ! uforce = uforce_others + uforce_boundary

  REAL :: vforce(nx,ny,nz)     ! forcing terms in v-momentum equ (kg/(m*s)**2)
                               ! vforce = vforce_others + vforce_boundary

  REAL :: wforce(nx,ny,nz)     ! forcing terms in w-momentum equ (kg/(m*s)**2)
                               ! wforce = wforce_others + wforce_boundary

  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model execution
  INCLUDE 'grid.inc'          ! Grid parameters
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 2, nx-1
        uforce(i,j,k) = uforce(i,j,k)                                   &
        - cbcdmp * 0.25 * (bcrlx(i-1,j)+bcrlx(i,j))  &
        * ( u(i,j,k) - ubar(i,j,k)) * (rhostr(i-1,j,k)+rhostr(i,j,k))
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 2, ny-1
      DO i = 1, nx-1
        vforce(i,j,k) = vforce(i,j,k)                                   &
        - cbcdmp * 0.25 * (bcrlx(i,j-1)+bcrlx(i,j))  &
        * ( v(i,j,k) - vbar(i,j,k)) * (rhostr(i,j-1,k)+rhostr(i,j,k))
      END DO
    END DO
  END DO

  DO k = 2, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        wforce(i,j,k) = wforce(i,j,k)    &
    - cbcdmp*bcrlx(i,j)*w(i,j,k)*0.5*(rhostr(i,j,k-1)+rhostr(i,j,k))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE brlxuvw_rbc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BRLXS_RBC                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE brlxs_rbc( nx,ny,nz,sprt,rhostr, ssrc, bcrlx )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the boundary relaxation and computational mixing for
!  a water/ice variable q in the boundary zone.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Xue Ming & Yuhe Liu
!  5/26/94
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    qprt     One of the water or ice variable
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!
!  INPUT/OUTPUT :
!
!    qsrc     Source terms in water/ice water equation (kg/(m*s)**2).
!             qvsrc = qvsrc_others + qvsrc_boundary
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc' ! Global constants that control model execution
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL :: sprt  (nx,ny,nz)       ! Water/ice variable
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: ssrc(nx,ny,nz)       ! Source term in q euation
  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k = 1, nz-1
    DO j = 1, ny
      DO i = 1, nx
        ssrc(i,j,k) = ssrc(i,j,k)  &
        - cbcdmp*bcrlx(i,j)*sprt(i,j,k)*rhostr(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE brlxs_rbc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE EXTBDTINI                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE extbdtini_rbc(nx,ny,nz, bcrlx)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in predicted variables from the first available external data
!  sets to calculate the linear time tendencies of the external data
!  set.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  8/10/94
!
!  MODIFICATION HISTORY:
!
!  08/30/1995 (Yuhe Liu)
!  Changed the initial boundary arrays, for restart run, from the
!  model arrays to the time interplated external boundary arrays.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!  nx         Number of grid points in the x-direction (east/west)
!  ny         Number of grid points in the y-direction (north/south)
!  nz         Number of grid points in the vertical
!
!  nx,ny,nz   Number of grid points in x, y, and z dir.
!
!  OUTPUT:
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
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations and COMMON blocks.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz       ! Number of grid points in x, y, and z dir.

  REAL :: bcrlx (nx,ny)     ! EXBC relaxation coefficients
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k, n, nq
  INTEGER :: iebc,iwbc,jnbc,jsbc,idist
!
!-----------------------------------------------------------------------
!
!  Declare the external boundary fields.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Make sure the external boundary fields match and fit the model
!  domain.
!
!-----------------------------------------------------------------------
!
  n = 2*ngbrz
!
!-----------------------------------------------------------------------
!
!  Fill the boundary relaxation varibles with zero.
!
!-----------------------------------------------------------------------
!
  DO j = 1, ny
    DO i = 1, nx
      bcrlx(i,j) = 0.
    END DO
  END DO

  iebc=1
!     iwbc=nx-1
!     jnbc=ny-1
  iwbc=nx-1+(nproc_x-1)*(nx-3)
  jnbc=ny-1+(nproc_y-1)*(ny-3)
  jsbc=1
  IF (mp_opt > 0) THEN
  END IF

  DO j = 1, ny-1
    DO i = 1, nx-1
!         idist=min(i-iebc,iwbc-i,jnbc-j,j-jsbc)
      idist=MIN(i+(loc_x-1)*(nx-3)-iebc,                                &
                iwbc-i-(loc_x-1)*(nx-3),                                &
                jnbc-j-(loc_y-1)*(ny-3),                                &
                j+(loc_y-1)*(ny-3)-jsbc)
      IF(idist < ngbrz ) bcrlx(i,j) = 1./(1.+(FLOAT(idist)/brlxhw)**2)
    END DO
  END DO

  IF (myproc == 0) THEN
    WRITE (6,'(/a)') 'Boundary relaxation coefficients'
    WRITE (6,'(/a,a)')                                                  &
        ' j\\i    1       2       3       4       5       6       7',   &
        '       8       9'
    WRITE (6,'(i3,9f8.5)') (j,(bcrlx(i,j),i=1,9),j=1,ny)
  END IF


  RETURN
END SUBROUTINE extbdtini_rbc
