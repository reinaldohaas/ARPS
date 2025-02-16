!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INCRSAVE                   ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######              University of Oklahoma.                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE incrsave(nx,ny,nz,nxndg,nyndg,nzndg,                         &
           u,v,w,pprt,ptprt,qv,                                         &
           qscalar,                                                     &
           uincr,vincr,wincr,pincr,ptincr,qvincr,                       &
           qcincr,qrincr,qiincr,qsincr,qhincr,                          &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:  Save the original background field in the incrememt
!  arrays for later computiona of the total analysis increment.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  March 1998
!
!  MODIFICATION HISTORY:
!
!  07/10/2001 (K. Brewster)
!  Changed to accomodate dynamically allocated increment arrays.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER :: nx,ny,nz
  INTEGER :: nxndg,nyndg,nzndg

  REAL :: u(nx,ny,nz)
  REAL :: v(nx,ny,nz)
  REAL :: w(nx,ny,nz)
  REAL :: pprt(nx,ny,nz)
  REAL :: ptprt(nx,ny,nz)
  REAL :: qv(nx,ny,nz)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: uincr(nxndg,nyndg,nzndg)
  REAL :: vincr(nxndg,nyndg,nzndg)
  REAL :: wincr(nxndg,nyndg,nzndg)
  REAL :: pincr(nxndg,nyndg,nzndg)
  REAL :: ptincr(nxndg,nyndg,nzndg)
  REAL :: qvincr(nxndg,nyndg,nzndg)
  REAL :: qcincr(nxndg,nyndg,nzndg)
  REAL :: qrincr(nxndg,nyndg,nzndg)
  REAL :: qiincr(nxndg,nyndg,nzndg)
  REAL :: qsincr(nxndg,nyndg,nzndg)
  REAL :: qhincr(nxndg,nyndg,nzndg)
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,nq
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nzndg
    DO j=1,nyndg
      DO i=1,nxndg

        uincr(i,j,k) =    u(i,j,k)
        vincr(i,j,k) =    v(i,j,k)
        wincr(i,j,k) =    w(i,j,k)
        pincr(i,j,k) = pprt(i,j,k)
        ptincr(i,j,k)=ptprt(i,j,k)
        qvincr(i,j,k)=   qv(i,j,k)
        qcincr(i,j,k)= 0.0
        qrincr(i,j,k)= 0.0
        qiincr(i,j,k)= 0.0
        qsincr(i,j,k)= 0.0
        qhincr(i,j,k)= 0.0

        IF (P_QC > 0) qcincr(i,j,k) = qscalar(i,j,k,P_QC)
        IF (P_QR > 0) qrincr(i,j,k) = qscalar(i,j,k,P_QR)
        IF (P_QI > 0) qiincr(i,j,k) = qscalar(i,j,k,P_QI)
        IF (P_QS > 0) qsincr(i,j,k) = qscalar(i,j,k,P_QS)
        IF (P_QH > 0) qhincr(i,j,k) = qscalar(i,j,k,P_QH)

      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE incrsave
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INCRCALC                   ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######              University of Oklahoma.                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE incrcalc(nx,ny,nz,nxndg,nyndg,nzndg,                         &
           u,v,w,pprt,ptprt,qv,                                         &
           qscalar,                                                     &
           uincr,vincr,wincr,pincr,ptincr,qvincr,                       &
           qcincr,qrincr,qiincr,qsincr,qhincr,                          &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE: Calculate total analysis increment by subtracting
!  the analyzed fields from the original background fields which
!  were temporarily held in the incr arrays.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  March 1998
!
!  MODIFICATION HISTORY:
!
!  07/10/2001 (K. Brewster)
!  Changed to accomodate dynamically allocated increment arrays.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER :: nx,ny,nz
  INTEGER :: nxndg,nyndg,nzndg

  REAL :: u(nx,ny,nz)
  REAL :: v(nx,ny,nz)
  REAL :: w(nx,ny,nz)
  REAL :: pprt(nx,ny,nz)
  REAL :: ptprt(nx,ny,nz)
  REAL :: qv(nx,ny,nz)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: uincr (nxndg,nyndg,nzndg)
  REAL :: vincr (nxndg,nyndg,nzndg)
  REAL :: wincr (nxndg,nyndg,nzndg)
  REAL :: pincr (nxndg,nyndg,nzndg)
  REAL :: ptincr(nxndg,nyndg,nzndg)
  REAL :: qvincr(nxndg,nyndg,nzndg)
  REAL :: qcincr(nxndg,nyndg,nzndg)
  REAL :: qrincr(nxndg,nyndg,nzndg)
  REAL :: qiincr(nxndg,nyndg,nzndg)
  REAL :: qsincr(nxndg,nyndg,nzndg)
  REAL :: qhincr(nxndg,nyndg,nzndg)
  INTEGER :: istatus

!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,nq
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nzndg
    DO j=1,nyndg
      DO i=1,nxndg

        uincr(i,j,k) =    u(i,j,k)- uincr(i,j,k)
        vincr(i,j,k) =    v(i,j,k)- vincr(i,j,k)
        wincr(i,j,k) =    w(i,j,k)- wincr(i,j,k)
        pincr(i,j,k) = pprt(i,j,k)- pincr(i,j,k)
        ptincr(i,j,k)=ptprt(i,j,k)-ptincr(i,j,k)
        qvincr(i,j,k)=   qv(i,j,k)-qvincr(i,j,k)
        qcincr(i,j,k)= 0.0
        qrincr(i,j,k)= 0.0
        qiincr(i,j,k)= 0.0
        qsincr(i,j,k)= 0.0
        qhincr(i,j,k)= 0.0
        IF (P_QC > 0) qcincr(i,j,k)= qscalar(i,j,k,P_QC)-qcincr(i,j,k)
        IF (P_QR > 0) qrincr(i,j,k)= qscalar(i,j,k,P_QR)-qrincr(i,j,k)
        IF (P_QI > 0) qiincr(i,j,k)= qscalar(i,j,k,P_QI)-qiincr(i,j,k)
        IF (P_QS > 0) qsincr(i,j,k)= qscalar(i,j,k,P_QS)-qsincr(i,j,k)
        IF (P_QH > 0) qhincr(i,j,k)= qscalar(i,j,k,P_QH)-qhincr(i,j,k)

      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE incrcalc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INCRDUMP                   ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######              University of Oklahoma.                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE incrdump(nxndg,nyndg,nzndg,incrdmp,incrcompr,incrfnam,       &
           uincr,vincr,wincr,pincr,ptincr,qvincr,                       &
           qcincr,qrincr,qiincr,qsincr,qhincr,                          &
           uincdmp,vincdmp,wincdmp,                                     &
           pincdmp,ptincdmp,qvincdmp,                                   &
           qcincdmp,qrincdmp,qiincdmp,qsincdmp,qhincdmp,                &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write analysis increments shift vectors to a file for use in
!  Newtonian nudging assimilation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  March, 1998
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nxndg,nyndg,nzndg
  INTEGER :: incrdmp         ! Option to write ADAS analysis increments
                             ! to a file
                             ! = 0, Don't create increment file;
                             ! = 1, write increment file
                             !      (unformatted IEEE binary);
                             ! = 3, write HDF format increment file.

  INTEGER :: incrcompr       ! HDF4 compressed option

  CHARACTER (LEN=*) :: incrfnam
  CHARACTER (LEN=256) :: savename
  REAL :: uincr(nxndg,nyndg,nzndg)
  REAL :: vincr(nxndg,nyndg,nzndg)
  REAL :: wincr(nxndg,nyndg,nzndg)
  REAL :: pincr(nxndg,nyndg,nzndg)
  REAL :: ptincr(nxndg,nyndg,nzndg)
  REAL :: qvincr(nxndg,nyndg,nzndg)
  REAL :: qcincr(nxndg,nyndg,nzndg)
  REAL :: qrincr(nxndg,nyndg,nzndg)
  REAL :: qiincr(nxndg,nyndg,nzndg)
  REAL :: qsincr(nxndg,nyndg,nzndg)
  REAL :: qhincr(nxndg,nyndg,nzndg)
!
  INTEGER :: uincdmp,vincdmp,wincdmp,pincdmp,ptincdmp,qvincdmp,         &
          qcincdmp,qrincdmp,qiincdmp,qsincdmp,qhincdmp
!
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i4time
  INTEGER :: nchinc,ierr

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:)
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, allocatable :: hmax(:), hmin(:) ! Temporary array
  INTEGER :: sd_id, stat

  INTEGER :: i,j,k
  INTEGER :: nxlgndg, nylgndg

  LOGICAL :: JOIN_FIELD
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'mp.inc'
  INCLUDE 'bndry.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (incrdmp == 0) RETURN

  IF (mp_opt > 0 .AND. joindmp(FINDX_I) == 0) THEN
    CALL gtsplitfn(incrfnam,1,1,loc_x,loc_y,1,1,0,0,0,2,savename,istatus)
  ELSE
    savename = incrfnam
  END IF

  JOIN_FIELD = .FALSE.
  IF (mp_opt > 0 .AND. joindmp(FINDX_I) /= 0 ) THEN
    nxlgndg = (nxndg - 3) * nproc_x + 3
    nylgndg = (nyndg - 3) * nproc_y + 3
    JOIN_FIELD = .TRUE.
  ELSE
    nxlgndg = nxndg
    nylgndg = nyndg
  END IF

  IF (mp_opt > 0) THEN
    ALLOCATE(tem1(nxlgndg,nylgndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "incrdump:tem1")
  END IF

  IF (incrcompr > 3) THEN
    allocate (itmp(nxndg,nyndg,nzndg),stat=stat)
    IF (stat /= 0) THEN
      WRITE (6,*) "INCRDUMP: ERROR allocating itmp, returning"
      RETURN
    END IF
    allocate (hmax(nzndg),stat=stat)
    IF (stat /= 0) THEN
      WRITE (6,*) "INCRDUMP: ERROR allocating hmax, returning"
      RETURN
    END IF
    allocate (hmin(nzndg),stat=stat)
    IF (stat /= 0) THEN
      WRITE (6,*) "INCRDUMP: ERROR allocating hmin, returning"
      RETURN
    END IF
  END IF

!-----------------------------------------------------------------------
!
!  Write out in Fortran unformatted.
!
!-----------------------------------------------------------------------

  IF (incrdmp == 1) THEN

    IF (myproc == 0 .OR. (mp_opt > 0 .AND. joindmp(FINDX_I) == 0)) THEN

      CALL getunit(nchinc)

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(savename,'-F f77 -N ieee', ierr)

      OPEN(nchinc,FILE=TRIM(savename),ERR=950,                        &
       FORM='unformatted',STATUS='unknown')

      GO TO 900
      950 CONTINUE

      WRITE(6,'(/a,a,a/)')                                            &
          'INCRDUMP: Error opening analysis increment output file: ', &
          TRIM(savename),'. No file was written.'
      istatus = -1
      CALL arpsstop("Problem opening INCR file.",1)

      RETURN
      900 CONTINUE

      CALL ctim2abss( year,month,day,hour,minute,second, i4time)

      WRITE(nchinc) runname,nxndg,nyndg,nzndg,i4time,                   &
                year,month,day,hour,minute,second

      WRITE(nchinc) mapproj,trulat1,trulat2,trulon,                     &
                sclfct,ctrlat,ctrlon

      WRITE(nchinc) uincdmp,vincdmp,wincdmp,pincdmp,ptincdmp,qvincdmp,  &
                    qcincdmp,qrincdmp,qiincdmp,qsincdmp,qhincdmp
    END IF

    IF(uincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(uincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'uinc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'uinc    '
        WRITE(nchinc) uincr
      END IF
    END IF

    IF(vincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(vincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'vinc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'vinc    '
        WRITE(nchinc) vincr
      END IF
    END IF

    IF(wincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(wincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'winc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'winc    '
        WRITE(nchinc) wincr
      END IF
    END IF

    IF(pincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(pincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'pinc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'pinc    '
        WRITE(nchinc) pincr
      END IF
    END IF

    IF(ptincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(ptincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'ptinc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'ptinc    '
        WRITE(nchinc) ptincr
      END IF
    END IF

    IF(qvincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qvincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'qvinc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'qvinc    '
        WRITE(nchinc) qvincr
      END IF
    END IF

    IF(qcincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qcincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'qcinc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'qcinc    '
        WRITE(nchinc) qcincr
      END IF
    END IF

    IF(qrincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qrincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'qrinc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'qrinc    '
        WRITE(nchinc) qrincr
      END IF
    END IF

    IF(qiincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qiincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'qiinc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'qiinc    '
        WRITE(nchinc) qiincr
      END IF
    END IF

    IF(qsincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qsincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'qsinc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'qsinc    '
        WRITE(nchinc) qsincr
      END IF
    END IF

    IF(qhincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qhincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0) THEN
          WRITE(nchinc) 'qhinc    '
          WRITE(nchinc) tem1
        ENDIF
      ELSE
        WRITE(nchinc) 'qhinc    '
        WRITE(nchinc) qhincr
      END IF
    END IF

!
    IF (myproc == 0 .OR. (mp_opt > 0 .AND. joindmp(FINDX_I) == 0)) THEN
      CLOSE(nchinc)
      CALL retunit(nchinc)
    END IF

!-----------------------------------------------------------------------
!
!  Write out in HDF4.
!
!-----------------------------------------------------------------------

  ELSE IF (incrdmp == 3) THEN


    IF (myproc == 0 .OR. (mp_opt > 0 .AND. joindmp(FINDX_I) == 0)) THEN

      CALL hdfopen(TRIM(savename), 2, sd_id)
      IF (sd_id < 0) THEN
        WRITE(6,'(/a,a,a/)')                                            &
          'INCRDUMP: Error opening analysis increment output file: ',   &
          TRIM(savename),'. No file was written.'
        istatus = -1
        CALL arpsstop("Problem opening INCR file.",1)
        RETURN
      END IF

      CALL ctim2abss( year,month,day,hour,minute,second, i4time)
      CALL hdfwrti(sd_id, 'i4time', i4time, stat)
      CALL hdfwrti(sd_id, 'day', day, stat)
      CALL hdfwrti(sd_id, 'year', year, stat)
      CALL hdfwrti(sd_id, 'month', month, stat)
      CALL hdfwrti(sd_id, 'hour', hour, stat)
      CALL hdfwrti(sd_id, 'minute', minute, stat)
      CALL hdfwrti(sd_id, 'second', second, stat)

      CALL hdfwrti(sd_id, 'nx', nxlgndg, stat)
      CALL hdfwrti(sd_id, 'ny', nylgndg, stat)
      CALL hdfwrti(sd_id, 'nz', nzndg, stat)
      CALL hdfwrtr(sd_id, 'dx', dx, stat)
      CALL hdfwrtr(sd_id, 'dy', dy, stat)
      CALL hdfwrtr(sd_id, 'dz', dz, stat)
      CALL hdfwrtr(sd_id, 'dzmin', dzmin, stat)
      CALL hdfwrti(sd_id, 'strhopt', strhopt, stat)
      CALL hdfwrtr(sd_id, 'zrefsfc', zrefsfc, stat)
      CALL hdfwrtr(sd_id, 'dlayer1', dlayer1, stat)
      CALL hdfwrtr(sd_id, 'dlayer2', dlayer2, stat)
      CALL hdfwrtr(sd_id, 'zflat', zflat, stat)
      CALL hdfwrtr(sd_id, 'strhtune', strhtune, stat)
      CALL hdfwrti(sd_id, 'mapproj', mapproj, stat)
      CALL hdfwrtr(sd_id, 'trulat1', trulat1, stat)
      CALL hdfwrtr(sd_id, 'trulat2', trulat2, stat)
      CALL hdfwrtr(sd_id, 'trulon', trulon, stat)
      CALL hdfwrtr(sd_id, 'sclfct', sclfct, stat)
      CALL hdfwrtr(sd_id, 'ctrlat', ctrlat, stat)
      CALL hdfwrtr(sd_id, 'ctrlon', ctrlon, stat)
    END IF

    IF(uincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(uincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)                                                &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
                  'uincr','u-velocity increment','m/s',                 &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(uincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,    &
                  'uincr','u-velocity increment','m/s',                 &
                   itmp,hmax,hmin)
      END IF
    END IF

    IF(vincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(vincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)                                                &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
                  'vincr','v-velocity increment','m/s',                 &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(vincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,    &
                  'vincr','v-velocity increment','m/s',                 &
                   itmp,hmax,hmin)
      END IF
    END IF

    IF(wincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(wincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)   &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
                  'wincr','w-velocity increment','m/s',                 &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(wincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,    &
                  'wincr','w-velocity increment','m/s',                 &
                   itmp,hmax,hmin)
      END IF
    END IF

    IF(pincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(pincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)                                                &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
                  'pincr','Pressure increment','Pascal',                &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(pincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,    &
                  'pincr','Pressure increment','Pascal',                &
                   itmp,hmax,hmin)
      END IF
    END IF

    IF(ptincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(ptincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)                                                &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
                  'ptincr','Potential temperature increment','K',       &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(ptincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,   &
                  'ptincr','Potential temperature increment','K',       &
                   itmp,hmax,hmin)
      END IF
    END IF

    IF(qvincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qvincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)                                                &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
          'qvincr','Water vapor specific humidity increment','kg/kg',   &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(qvincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,   &
          'qvincr','Water vapor specific humidity increment','kg/kg',   &
                   itmp,hmax,hmin)
      END IF
    END IF

    IF(qcincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qcincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)                                                &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
             'qcincr','Cloud water mixing ratio increment','kg/kg',     &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(qcincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,   &
             'qcincr','Cloud water mixing ratio increment','kg/kg',     &
                   itmp,hmax,hmin)
      END IF
    END IF

    IF(qrincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qrincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)                                                &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
              'qrincr','Rain water mixing ratio increment','kg/kg',     &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(qrincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,   &
              'qrincr','Rain water mixing ratio increment','kg/kg',     &
                   itmp,hmax,hmin)
      END IF
    END IF

    IF(qiincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qiincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)                                                &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
              'qiincr','Cloud ice mixing ratio increment','kg/kg',      &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(qiincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,   &
              'qiincr','Cloud ice mixing ratio increment','kg/kg',      &
                   itmp,hmax,hmin)
      END IF
    END IF

    IF(qsincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qsincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)                                                &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
                  'qsincr','Snow mixing ratio increment','kg/kg',       &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(qsincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,   &
                  'qsincr','Snow mixing ratio increment','kg/kg',       &
                   itmp,hmax,hmin)
      END IF

    END IF
    IF(qhincdmp > 0) THEN
      IF (JOIN_FIELD) THEN
        CALL mpimerge3d(qhincr,nxndg,nyndg,nzndg,tem1)
        IF (myproc == 0)                                                &
        CALL hdfwrt3d(tem1,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,     &
                  'qhincr','Hail mixing ratio increment','kg/kg',       &
                   itmp,hmax,hmin)
      ELSE
        CALL hdfwrt3d(qhincr,nxlgndg,nylgndg,nzndg,sd_id,1,incrcompr,   &
                  'qhincr','Hail mixing ratio increment','kg/kg',       &
                   itmp,hmax,hmin)
      END IF
    END IF

    IF (myproc == 0 .OR. (mp_opt > 0 .AND. joindmp(FINDX_I) == 0)) THEN
      CALL hdfclose(sd_id,stat)
      deallocate (itmp,stat=stat)
      deallocate (hmax,stat=stat)
      deallocate (hmin,stat=stat)
    END IF

  ELSE

    ! alternate dump format ...
    IF (myproc == 0)                                                    &
      WRITE(6,*) 'The supported increment data dump format are ',       &
               'binary (incrdmp=1) and HDF4 (incrdmp = 3).'
    CALL arpsstop('Increment data dump format is not supported.',1)

  END IF

  IF (myproc == 0)  WRITE(6,'(/a,a,a/)')                                &
      'INCRDUMP: Successfully wrote analysis incr file: ',              &
      TRIM(savename),'.'
  istatus = 0

  IF (ALLOCATED(tem1))  DEALLOCATE(tem1)

  RETURN
END SUBROUTINE incrdump
