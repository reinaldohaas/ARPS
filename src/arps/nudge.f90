!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ININUDGE                   ######
!######                                                      ######
!######                   Developed by                       ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ininudge(nxndg,nyndg,nzndg,                                 &
                    uincr,vincr,wincr,pincr,ptincr,qvincr,             &
                    qcincr,qrincr,qiincr,qsincr,qhincr,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize analysis increments for use in continuous
!  nudging adjustment process.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  March, 1998
!
!  MODIFICATION HISTORY:
!
!  07/10/2001 (K. Brewster)
!  Added increment arrays to argument list and removed
!  initialization to zero before reading (now done prior to call).
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

  REAL :: uincr(nxndg,nyndg,nzndg)      ! Analysis increment for u
  REAL :: vincr(nxndg,nyndg,nzndg)      ! Analysis increment for v
  REAL :: wincr(nxndg,nyndg,nzndg)      ! Analysis increment for w
  REAL :: pincr(nxndg,nyndg,nzndg)      ! Analysis increment for p
  REAL :: ptincr(nxndg,nyndg,nzndg)     ! Analysis increment for pt
  REAL :: qvincr(nxndg,nyndg,nzndg)     ! Analysis increment for qv
  REAL :: qcincr(nxndg,nyndg,nzndg)     ! Analysis increment for qc
  REAL :: qrincr(nxndg,nyndg,nzndg)     ! Analysis increment for qr
  REAL :: qiincr(nxndg,nyndg,nzndg)     ! Analysis increment for qi
  REAL :: qsincr(nxndg,nyndg,nzndg)     ! Analysis increment for qs
  REAL :: qhincr(nxndg,nyndg,nzndg)     ! Analysis increment for qh

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ncorx
  INTEGER :: i,j,k
!  real ndtime
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'nudging.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Initializations
!
!-----------------------------------------------------------------------
!
  IF ( nudgopt <= 0 ) THEN
    RETURN
  END IF

  ndtime=ndstop-ndstart
  IF (myproc == 0) WRITE(6,'(a,f10.2,a)')                               &
      ' Application of nudging adjustment lasts ',ndtime,' secs'
!
!-----------------------------------------------------------------------
!
!  Read analysis increment file
!
!-----------------------------------------------------------------------
!
  IF(mp_opt > 0 .AND. readsplit(FINDX_I) > 0) THEN
    CALL incrreadsplit(nxndg,nyndg,nzndg,incrfnam,                      &
                uincr,vincr,wincr,pincr,ptincr,qvincr,                  &
                qcincr,qrincr,qiincr,qsincr,qhincr,                     &
                istatus)
  ELSE
    CALL incrread(nxndg,nyndg,nzndg,incrfnam,                           &
                uincr,vincr,wincr,pincr,ptincr,qvincr,                  &
                qcincr,qrincr,qiincr,qsincr,qhincr,                     &
                istatus)
  END IF
!
!-----------------------------------------------------------------------
!
!  Compute the fixed time scale factor.
!
!-----------------------------------------------------------------------
!
  IF(ndintvl > 0.) THEN
    ncorx=nint(ndtime/ndintvl)
    IF (myproc == 0) WRITE(6,'(a,i5,a,/a,f10.2)')                       &
            ' Nudging applied in ',ncorx,' steps',                      &
            ' ndintvl adjusted for dtbig =',ndintvl
  ELSE
    WRITE(6,'(a,/a,f10.2)')                                             &
            ' Try again using new ndintvl',                             &
            ' ndintvl adjusted for dtbig =',ndintvl
    WRITE(6,'(a)') ' STOPPING in ININUDGE'
    CALL arpsstop('arpsstop called from ININUDGE improper nudging       &
                 & interval selected.',1)
  END IF

  RETURN
END SUBROUTINE ininudge
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INCRREAD                   ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######              University of Oklahoma.                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE incrread(nxndg,nyndg,nzndg,incrfnam,                         &
           uincr,vincr,wincr,pincr,ptincr,qvincr,                       &
           qcincr,qrincr,qiincr,qsincr,qhincr,                          &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read analysis increments from a file for use in continuous
!  adjustment process.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  March 1998
!
!  MODIFICATION HISTORY:
!
!  07/10/2001 (K. Brewster)
!  Added increment arrays to argument list rather than obtaining
!  them through common in nudging.inc
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
  CHARACTER (LEN=256) :: incrfnam

  REAL :: uincr(nxndg,nyndg,nzndg)      ! Analysis increment for u
  REAL :: vincr(nxndg,nyndg,nzndg)      ! Analysis increment for v
  REAL :: wincr(nxndg,nyndg,nzndg)      ! Analysis increment for w
  REAL :: pincr(nxndg,nyndg,nzndg)      ! Analysis increment for p
  REAL :: ptincr(nxndg,nyndg,nzndg)     ! Analysis increment for pt
  REAL :: qvincr(nxndg,nyndg,nzndg)     ! Analysis increment for qv
  REAL :: qcincr(nxndg,nyndg,nzndg)     ! Analysis increment for qc
  REAL :: qrincr(nxndg,nyndg,nzndg)     ! Analysis increment for qr
  REAL :: qiincr(nxndg,nyndg,nzndg)     ! Analysis increment for qi
  REAL :: qsincr(nxndg,nyndg,nzndg)     ! Analysis increment for qs
  REAL :: qhincr(nxndg,nyndg,nzndg)     ! Analysis increment for qh

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80) :: runnamin
  CHARACTER (LEN=8) :: varin
  INTEGER :: iyr,imon,idy,ihr,imin,isec
  INTEGER :: nxin,nyin,nzin,i4timein
  INTEGER :: maprojin
  REAL :: trlat1in,trlat2in,trlonin
  REAL :: sclfctin,ctrlatin,ctrlonin
  INTEGER :: ustor,vstor,wstor,pstor,ptstor,qvstor,                     &
          qcstor,qrstor,qistor,qsstor,qhstor

  INTEGER :: nchinc,ierr

  INTEGER :: ireturn
  INTEGER :: strhoptin
  REAL :: dxin,dyin,dzin,dzminin,zrefsfcin,dlayer1in,                   &
          dlayer2in,zflatin,strhtunein

  INTEGER(2), allocatable :: itmp(:,:,:) ! Temporary array
  REAL, allocatable :: hmax(:), hmin(:) ! Temporary array
  INTEGER :: istat, sd_id

  CHARACTER (LEN=256) :: savename
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (incrfmt == 3) THEN
    ALLOCATE (itmp(nxndg,nyndg,nzndg),stat=istat)
    CALL check_alloc_status(istat, "INCREAD:itmp")

    ALLOCATE (hmax(nzndg),stat=istat)
    CALL check_alloc_status(istat, "INCREAD:hmax")

    ALLOCATE (hmin(nzndg),stat=istat)
    CALL check_alloc_status(istat, "INCREAD:hmin")
  END IF

  !wdt kwthomas update
  IF (mp_opt > 0) THEN
    savename(1:256) = incrfnam(1:256)
    CALL gtsplitfn(savename,1,1,loc_x,loc_y,1,1,0,0,1,lvldbg,incrfnam,istat)
  END IF
!
!-----------------------------------------------------------------------
!
!  Get unit number and open file
!
!-----------------------------------------------------------------------
!
  IF (incrfmt == 1) THEN

!-----------------------------------------------------------------------
!
!  Fortran unformatted dump.
!
!-----------------------------------------------------------------------

    CALL getunit(nchinc)

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(incrfnam,'-F f77 -N ieee', ierr)

    OPEN(nchinc,FILE=trim(incrfnam),ERR=950,                            &
        FORM='unformatted',STATUS='old')

    READ(nchinc,ERR=950) runnamin,nxin,nyin,nzin,i4timein,              &
                iyr,imon,idy,ihr,imin,isec

    READ(nchinc,ERR=950) maprojin,trlat1in,trlat2in,trlonin,            &
                sclfctin,ctrlatin,ctrlonin

    READ(nchinc,ERR=950) ustor,vstor,wstor,pstor,ptstor,qvstor,         &
               qcstor,qrstor,qistor,qsstor,qhstor

    IF(ustor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into u-increment field'
      READ(nchinc,ERR=950) uincr
    END IF
    IF(vstor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into v-increment field'
      READ(nchinc,ERR=950) vincr
    END IF
    IF(wstor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into w-increment field'
      READ(nchinc,ERR=950) wincr
    END IF
    IF(pstor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into p-increment field'
      READ(nchinc,ERR=950) pincr
    END IF
    IF(ptstor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into pt-increment field'
      READ(nchinc,ERR=950) ptincr
    END IF
    IF(qvstor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into qv-increment field'
      READ(nchinc,ERR=950) qvincr
    END IF
    IF(qcstor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into qc-increment field'
      READ(nchinc,ERR=950) qcincr
    END IF
    IF(qrstor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into qr-increment field'
      READ(nchinc,ERR=950) qrincr
    END IF
    IF(qistor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into qi-increment field'
      READ(nchinc,ERR=950) qiincr
    END IF
    IF(qsstor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into qs-increment field'
      READ(nchinc,ERR=950) qsincr
    END IF
    IF(qhstor > 0) THEN
      READ(nchinc,ERR=950) varin
      WRITE(6,'(a,a,a)') 'Reading ',varin,' into qh-increment field'
      READ(nchinc,ERR=950) qhincr
    END IF

    WRITE(6,'(/a,a/)')                                                  &
        ' Successfully read analysis incr file: ',incrfnam
    istatus=0
    !wdt kwthomas update
    IF (mp_opt > 0) incrfnam(1:256) = savename(1:256)
    RETURN
  ELSE IF (incrfmt == 3) THEN

!-----------------------------------------------------------------------
!
!  HDF4 format.
!
!-----------------------------------------------------------------------

    CALL hdfopen(trim(incrfnam), 1, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "INCRREAD: ERROR opening ",                           &
                 trim(incrfnam)," for reading."
      istatus = 1
      GO TO 900
    END IF

    CALL hdfrdi(sd_id,"i4time",i4timein,istat)
    CALL hdfrdi(sd_id,"month",imon,istat)
    CALL hdfrdi(sd_id,"day",idy,istat)
    CALL hdfrdi(sd_id,"year",iyr,istat)
    CALL hdfrdi(sd_id,"hour",ihr,istat)
    CALL hdfrdi(sd_id,"minute",imin,istat)
    CALL hdfrdi(sd_id,"second",isec,istat)

    CALL hdfrdi(sd_id,"nx",nxin,istat)
    CALL hdfrdi(sd_id,"ny",nyin,istat)
    CALL hdfrdi(sd_id,"nz",nzin,istat)
    CALL hdfrdr(sd_id,"dx",dxin,istat)
    CALL hdfrdr(sd_id,"dy",dyin,istat)
    CALL hdfrdr(sd_id,"dz",dzin,istat)
    CALL hdfrdr(sd_id,"dzmin",dzminin,istat)
    CALL hdfrdi(sd_id,"strhopt",strhoptin,istat)
    CALL hdfrdr(sd_id,"zrefsfc",zrefsfcin,istat)
    CALL hdfrdr(sd_id,"dlayer1",dlayer1in,istat)
    CALL hdfrdr(sd_id,"dlayer2",dlayer2in,istat)
    CALL hdfrdr(sd_id,"zflat",zflatin,istat)
    CALL hdfrdr(sd_id,"strhtune",strhtunein,istat)
    CALL hdfrdi(sd_id,"mapproj",maprojin,istat)
    CALL hdfrdr(sd_id,"trulat1",trlat1in,istat)
    CALL hdfrdr(sd_id,"trulat2",trlat2in,istat)
    CALL hdfrdr(sd_id,"trulon",trlonin,istat)
    CALL hdfrdr(sd_id,"sclfct",sclfctin,istat)
    CALL hdfrdr(sd_id,"ctrlat",ctrlatin,istat)
    CALL hdfrdr(sd_id,"ctrlon",ctrlonin,istat)

    CALL checkgrid3d(nxndg,nyndg,nzndg,nxin,nyin,nzin,                  &
        dx,dy,dz,dzmin,ctrlat,ctrlon,                                   &
        strhopt,zrefsfc,dlayer1,dlayer2,zflat,strhtune,                 &
        mapproj,trulat1,trulat2,trulon,sclfct,                          &
        dxin,dyin,dzin,dzminin,ctrlatin,ctrlonin,                       &
        strhoptin,zrefsfcin,dlayer1in,dlayer2in,zflatin,strhtunein,     &
        maprojin,trlat1in,trlat2in,trlonin,sclfctin,ireturn)

    IF (ireturn /= 0) THEN
      WRITE (6,*) "INCRREAD: ERROR, grid parameter mismatch"
      istatus = 1
      GO TO 900
    END IF

    CALL hdfrdi(sd_id,"i4time",i4timein,istat)
    CALL hdfrdi(sd_id,"iyr",iyr,istat)
    CALL hdfrdi(sd_id,"imon",imon,istat)
    CALL hdfrdi(sd_id,"idy",idy,istat)
    CALL hdfrdi(sd_id,"ihr",ihr,istat)
    CALL hdfrdi(sd_id,"imin",imin,istat)
    CALL hdfrdi(sd_id,"isec",isec,istat)

    CALL hdfrd3d(sd_id,"uincr",nxndg,nyndg,nzndg,uincr,                 &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950
    CALL hdfrd3d(sd_id,"vincr",nxndg,nyndg,nzndg,vincr,                 &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950
    CALL hdfrd3d(sd_id,"wincr",nxndg,nyndg,nzndg,wincr,                 &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950
    CALL hdfrd3d(sd_id,"pincr",nxndg,nyndg,nzndg,pincr,                 &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950
    CALL hdfrd3d(sd_id,"ptincr",nxndg,nyndg,nzndg,ptincr,               &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950
    CALL hdfrd3d(sd_id,"qvincr",nxndg,nyndg,nzndg,qvincr,               &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950
    CALL hdfrd3d(sd_id,"qcincr",nxndg,nyndg,nzndg,qcincr,               &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950
    CALL hdfrd3d(sd_id,"qrincr",nxndg,nyndg,nzndg,qrincr,               &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950
    CALL hdfrd3d(sd_id,"qiincr",nxndg,nyndg,nzndg,qiincr,               &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950
    CALL hdfrd3d(sd_id,"qsincr",nxndg,nyndg,nzndg,qsincr,               &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950
    CALL hdfrd3d(sd_id,"qhincr",nxndg,nyndg,nzndg,qhincr,               &
                 istat,itmp,hmax,hmin)
    IF (istat > 1) GO TO 950

    istatus = 0

  ELSE

    ! alternate dump format ...
    WRITE(6,*) 'The supported increment data format are ',                &
               'binary (incrfmt=1) and HDF4 no compressed (incrfmt=3).'
    CALL arpsstop('Increment data format is not supported.',1)

  END IF

  GO TO 900

  950   CONTINUE
  WRITE(6,'(/a,a/)')                                                    &
      'INCRREAD: Error reading analysis incr output file: ',            &
      trim(incrfnam)
  istatus = 1
  WRITE(6,*) "INCCREAD: calling arpsstop"
  CALL arpsstop('arpsstop called from INCREAD error reading incr        &
               & output file.',1)

  900   CONTINUE

  IF (mp_opt > 0) incrfnam(1:256) = savename(1:256)

  IF (incrfmt == 1) THEN
    CLOSE(nchinc)
  ELSE
    CALL hdfclose(sd_id,istat)

    DEALLOCATE (itmp,stat=istat)
    DEALLOCATE (hmax,stat=istat)
    DEALLOCATE (hmin,stat=istat)
  END IF

  RETURN

END SUBROUTINE incrread

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INCRREADSPLIT              ######
!######                                                      ######
!######    Center for Analysis and Prediction of Storms      ######
!######              University of Oklahoma.                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE incrreadsplit(nxndg,nyndg,nzndg,incrfnam,                    &
           uincr,vincr,wincr,pincr,ptincr,qvincr,                       &
           qcincr,qrincr,qiincr,qsincr,qhincr,                          &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read and split analysis increments from a file for use in continuous
!  adjustment process. This subroutine is for mpi runs and it is base
!  on subroutine incrread.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  November 2002
!
!  MODIFICATION HISTORY:
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
  CHARACTER (LEN=256) :: incrfnam

  REAL :: uincr(nxndg,nyndg,nzndg)      ! Analysis increment for u
  REAL :: vincr(nxndg,nyndg,nzndg)      ! Analysis increment for v
  REAL :: wincr(nxndg,nyndg,nzndg)      ! Analysis increment for w
  REAL :: pincr(nxndg,nyndg,nzndg)      ! Analysis increment for p
  REAL :: ptincr(nxndg,nyndg,nzndg)     ! Analysis increment for pt
  REAL :: qvincr(nxndg,nyndg,nzndg)     ! Analysis increment for qv
  REAL :: qcincr(nxndg,nyndg,nzndg)     ! Analysis increment for qc
  REAL :: qrincr(nxndg,nyndg,nzndg)     ! Analysis increment for qr
  REAL :: qiincr(nxndg,nyndg,nzndg)     ! Analysis increment for qi
  REAL :: qsincr(nxndg,nyndg,nzndg)     ! Analysis increment for qs
  REAL :: qhincr(nxndg,nyndg,nzndg)     ! Analysis increment for qh

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80) :: runnamin
  CHARACTER (LEN=8) :: varin
  INTEGER :: iyr,imon,idy,ihr,imin,isec
  INTEGER :: nxin,nyin,nzin,i4timein
  INTEGER :: maprojin
  REAL :: trlat1in,trlat2in,trlonin
  REAL :: sclfctin,ctrlatin,ctrlonin
  INTEGER :: ustor,vstor,wstor,pstor,ptstor,qvstor,                     &
          qcstor,qrstor,qistor,qsstor,qhstor

  INTEGER :: nchinc,ierr

  INTEGER :: ireturn
  INTEGER :: strhoptin
  REAL :: dxin,dyin,dzin,dzminin,zrefsfcin,dlayer1in,                   &
          dlayer2in,zflatin,strhtunein

  INTEGER(2), allocatable :: itmp(:,:,:) ! Temporary array
  REAL, allocatable :: hmax(:), hmin(:) ! Temporary array
  INTEGER :: istat, sd_id

  INTEGER :: nxndglg, nyndglg
  REAL, ALLOCATABLE :: var3din(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxndglg = (nxndg-3)*nproc_x+3
  nyndglg = (nyndg-3)*nproc_y+3

  ALLOCATE(var3din(nxndglg,nyndglg,nzndg),stat=istat)
  CALL check_alloc_status(istat, "INCREADSPLIT:var3din")

  IF (incrfmt == 3) THEN
    ALLOCATE (itmp(nxndglg,nyndglg,nzndg),stat=istat)
    CALL check_alloc_status(istat, "INCREADSPLIT:itmp")

    ALLOCATE (hmax(nzndg),stat=istat)
    CALL check_alloc_status(istat, "INCREADSPLIT:itmp")

    ALLOCATE (hmin(nzndg),stat=istat)
    CALL check_alloc_status(istat, "INCREADSPLIT:hmin")
  END IF

!
!-----------------------------------------------------------------------
!
!  Get unit number and open file
!
!-----------------------------------------------------------------------
!
  IF (incrfmt == 1) THEN

!-----------------------------------------------------------------------
!
!  Fortran unformatted dump.
!
!-----------------------------------------------------------------------

    IF (myproc == 0) THEN
      CALL getunit(nchinc)

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(incrfnam,'-F f77 -N ieee', ierr)

      OPEN(nchinc,FILE=trim(incrfnam),ERR=950,                          &
          FORM='unformatted',STATUS='old')

      READ(nchinc,ERR=950) runnamin,nxin,nyin,nzin,i4timein,            &
                  iyr,imon,idy,ihr,imin,isec

      READ(nchinc,ERR=950) maprojin,trlat1in,trlat2in,trlonin,          &
                  sclfctin,ctrlatin,ctrlonin

      READ(nchinc,ERR=950) ustor,vstor,wstor,pstor,ptstor,qvstor,       &
                 qcstor,qrstor,qistor,qsstor,qhstor
    END IF

    CALL mpupdatei(nxin,1)          ! need to check nxin, and nyin
    CALL mpupdatei(nyin,1)          ! to do later -- WYH.
    CALL mpupdatei(nzin,1)

    CALL mpupdatec(runnamin,80)
    CALL mpupdatei(i4timein,1)
    CALL mpupdatei(iyr,1)
    CALL mpupdatei(imon,1)
    CALL mpupdatei(idy,1)
    CALL mpupdatei(ihr,1)
    CALL mpupdatei(imin,1)
    CALL mpupdatei(isec,1)

    CALL mpupdatei(maprojin,1)
    CALL mpupdater(trlat1in,1)
    CALL mpupdater(trlat2in,1)
    CALL mpupdater(trlonin,1)
    CALL mpupdater(sclfctin,1)
    CALL mpupdater(ctrlatin,1)
    CALL mpupdater(ctrlonin,1)

    CALL mpupdater(ustor,1)
    CALL mpupdater(vstor,1)
    CALL mpupdater(wstor,1)
    CALL mpupdater(pstor,1)
    CALL mpupdater(ptstor,1)
    CALL mpupdater(qvstor,1)
    CALL mpupdater(qcstor,1)
    CALL mpupdater(qrstor,1)
    CALL mpupdater(qistor,1)
    CALL mpupdater(qsstor,1)
    CALL mpupdater(qhstor,1)

    IF(ustor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into u-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,uincr)
    END IF

    IF(vstor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into v-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,vincr)
    END IF

    IF(wstor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into w-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,wincr)
    END IF

    IF(pstor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into p-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,pincr)
    END IF

    IF(ptstor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into pt-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,ptincr)
    END IF

    IF(qvstor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into qv-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qvincr)
    END IF

    IF(qcstor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into qc-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qcincr)
    END IF

    IF(qrstor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into qr-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qrincr)
    END IF

    IF(qistor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into qi-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qiincr)
    END IF

    IF(qsstor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into qs-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qsincr)
    END IF

    IF(qhstor > 0) THEN
      IF(myproc ==0) THEN
        READ(nchinc,ERR=950) varin
        WRITE(6,'(a,a,a)') 'Reading ',varin,' into qh-increment field'
        READ(nchinc,ERR=950) var3din
      END IF
      CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qhincr)
    END IF

  ELSE IF (incrfmt == 3) THEN

!-----------------------------------------------------------------------
!
!  HDF4 format.
!
!-----------------------------------------------------------------------

    IF(myproc ==0) THEN
      CALL hdfopen(trim(incrfnam), 1, sd_id)
      IF (sd_id < 0) THEN
        WRITE (6,*) "INCRREADSPLIT: ERROR opening ",                    &
                   trim(incrfnam)," for reading."
        istatus = 1
        GO TO 900
      END IF

      CALL hdfrdi(sd_id,"i4time",i4timein,istat)
      CALL hdfrdi(sd_id,"month",imon,istat)
      CALL hdfrdi(sd_id,"day",idy,istat)
      CALL hdfrdi(sd_id,"year",iyr,istat)
      CALL hdfrdi(sd_id,"hour",ihr,istat)
      CALL hdfrdi(sd_id,"minute",imin,istat)
      CALL hdfrdi(sd_id,"second",isec,istat)

      CALL hdfrdi(sd_id,"nx",nxin,istat)
      CALL hdfrdi(sd_id,"ny",nyin,istat)
      CALL hdfrdi(sd_id,"nz",nzin,istat)
      CALL hdfrdr(sd_id,"dx",dxin,istat)
      CALL hdfrdr(sd_id,"dy",dyin,istat)
      CALL hdfrdr(sd_id,"dz",dzin,istat)
      CALL hdfrdr(sd_id,"dzmin",dzminin,istat)
      CALL hdfrdi(sd_id,"strhopt",strhoptin,istat)
      CALL hdfrdr(sd_id,"zrefsfc",zrefsfcin,istat)
      CALL hdfrdr(sd_id,"dlayer1",dlayer1in,istat)
      CALL hdfrdr(sd_id,"dlayer2",dlayer2in,istat)
      CALL hdfrdr(sd_id,"zflat",zflatin,istat)
      CALL hdfrdr(sd_id,"strhtune",strhtunein,istat)
      CALL hdfrdi(sd_id,"mapproj",maprojin,istat)
      CALL hdfrdr(sd_id,"trulat1",trlat1in,istat)
      CALL hdfrdr(sd_id,"trulat2",trlat2in,istat)
      CALL hdfrdr(sd_id,"trulon",trlonin,istat)
      CALL hdfrdr(sd_id,"sclfct",sclfctin,istat)
      CALL hdfrdr(sd_id,"ctrlat",ctrlatin,istat)
      CALL hdfrdr(sd_id,"ctrlon",ctrlonin,istat)
    END IF

    CALL mpupdatei(i4timein,1)
    CALL mpupdatei(iyr,1)
    CALL mpupdatei(imon,1)
    CALL mpupdatei(idy,1)
    CALL mpupdatei(ihr,1)
    CALL mpupdatei(imin,1)
    CALL mpupdatei(isec,1)

    CALL mpupdatei(nxin,1)
    CALL mpupdatei(nyin,1)
    CALL mpupdatei(nzin,1)

    CALL mpupdater(dxin,1)
    CALL mpupdater(dyin,1)
    CALL mpupdater(dzin,1)
    CALL mpupdater(dzminin,1)
    CALL mpupdatei(strhoptin,1)
    CALL mpupdater(zrefsfcin,1)
    CALL mpupdater(dlayer1in,1)
    CALL mpupdater(dlayer2in,1)
    CALL mpupdater(zflatin,1)
    CALL mpupdater(strhtunein,1)

    CALL mpupdatei(maprojin,1)
    CALL mpupdater(trlat1in,1)
    CALL mpupdater(trlat2in,1)
    CALL mpupdater(trlonin,1)
    CALL mpupdater(sclfctin,1)
    CALL mpupdater(ctrlatin,1)
    CALL mpupdater(ctrlonin,1)

    CALL checkgrid3d(nxndglg,nyndglg,nzndg,nxin,nyin,nzin,              &
          dx,dy,dz,dzmin,ctrlat,ctrlon,                                 &
          strhopt,zrefsfc,dlayer1,dlayer2,zflat,strhtune,               &
          mapproj,trulat1,trulat2,trulon,sclfct,                        &
          dxin,dyin,dzin,dzminin,ctrlatin,ctrlonin,                     &
          strhoptin,zrefsfcin,dlayer1in,dlayer2in,zflatin,strhtunein,   &
          maprojin,trlat1in,trlat2in,trlonin,sclfctin,ireturn)

    IF (ireturn /= 0) THEN
      WRITE (6,*) "INCRREADSPLIT: ERROR, grid parameter mismatch"
      istatus = 1
      GO TO 900
    END IF

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"uincr",nxndglg,nyndglg,nzndg,var3din,        &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,uincr)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"vincr",nxndglg,nyndglg,nzndg,var3din,       &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,vincr)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"wincr",nxndglg,nyndglg,nzndg,var3din,       &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,wincr)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"pincr",nxndglg,nyndglg,nzndg,var3din,      &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,pincr)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"ptincr",nxndglg,nyndglg,nzndg,var3din,     &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,ptincr)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"qvincr",nxndglg,nyndglg,nzndg,var3din,    &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qvincr)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"qcincr",nxndglg,nyndglg,nzndg,var3din,    &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qcincr)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"qrincr",nxndglg,nyndglg,nzndg,var3din,   &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qrincr)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"qiincr",nxndglg,nyndglg,nzndg,var3din,   &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qiincr)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"qsincr",nxndglg,nyndglg,nzndg,var3din,   &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qsincr)

    IF(myproc == 0) THEN
      CALL hdfrd3d(sd_id,"qhincr",nxndglg,nyndglg,nzndg,var3din,   &
                   istat,itmp,hmax,hmin)
      IF (istat > 1) GO TO 950
    END IF
    CALL mpisplit3d(var3din,nxndg,nyndg,nzndg,qhincr)

  ELSE

    ! alternate data format ...
    WRITE(6,*) 'The supported Increment data format are ',           &
               'binary (incrfmt=1) and HDF4 no compressed (incrfmt = 3).'
    CALL arpsstop('Increment data format is not supported.',1)

  END IF

  IF(myproc ==0) WRITE(6,'(/a,a/)')                                  &
        ' Successfully read analysis incr file: ',incrfnam
  istatus=0

  GO TO 900

  950   CONTINUE
  WRITE(6,'(/a,a/)')                                                 &
      'INCRREADSPLIT: Error reading analysis incr output file: ',    &
      trim(incrfnam)
  istatus = 1

  WRITE(6,*) "INCRREADSPLIT: calling arpsstop"
  CALL arpsstop('arpsstop called from INCREAD error reading incr     &
               & output file.',1)

  900   CONTINUE

  IF(myproc == 0) THEN
    IF (incrfmt == 1) THEN
      CLOSE(nchinc)
    ELSE IF (incrfmt == 3) THEN
      CALL hdfclose(sd_id,istat)
    ELSE
      CALL arpsstop('Unknown formt in INCRREADSPLIT.',1)
    END IF
  END IF

  DEALLOCATE (var3din, stat=istat)
  IF (incrfmt == 3) DEALLOCATE (itmp,hmax,hmin,stat=istat)

  RETURN
END SUBROUTINE incrreadsplit

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE NUDGEALL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE nudgeall(nx,ny,nz,nxndg,nyndg,nzndg,                          &
           u,v,w,pprt,ptprt,qv,qscalar,                                  &
           uincr,vincr,wincr,pincr,ptincr,qvincr,                        &
           qcincr,qrincr,qiincr,qsincr,qhincr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  March, 1998
!
!  MODIFICATION HISTORY:
!
!  07/10/2001 (K. Brewster)
!  Added increment arrays to argument list rather than obtaining
!  them through common in nudging.inc
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nxndg    Number x grid points for IAU
!    nyndg    Number y grid points for IAU
!    nzndg    Number z grid points for IAU
!
!    u        x component of velocity at times tpast and tpresent (m/s)
!    v        y component of velocity at times tpast and tpresent (m/s)
!    w        Vertical component of Cartesian velocity at times
!             tpast and tpresent (m/s)
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and tpresent (Pascal)
!    qv       Water vapor specific humidity at times tpast and tpresent (kg/kg)
!    qc       Cloud water mixing ratio at times tpast and tpresent (kg/kg)
!    qr       Rainwater mixing ratio at times tpast and tpresent (kg/kg)
!    qi       Cloud ice mixing ratio at times tpast and tpresent (kg/kg)
!    qs       Snow mixing ratio at times tpast and tpresent (kg/kg)
!    qh       Hail mixing ratio at times tpast and tpresent (kg/kg)
!
!-----------------------------------------------------------------------

  IMPLICIT NONE             ! Force explicit declarations
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'nudging.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  INTEGER :: nxndg,nyndg,nzndg ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)  ! Perturbation potential temperature
                            ! from that of base state atmosphere (K)
  REAL :: pprt  (nx,ny,nz)  ! Perturbation pressure from that
                            ! of base state atmosphere (Pascal)
  REAL :: qv    (nx,ny,nz)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: uincr(nxndg,nyndg,nzndg)      ! Analysis increment for u
  REAL :: vincr(nxndg,nyndg,nzndg)      ! Analysis increment for v
  REAL :: wincr(nxndg,nyndg,nzndg)      ! Analysis increment for w
  REAL :: pincr(nxndg,nyndg,nzndg)      ! Analysis increment for p
  REAL :: ptincr(nxndg,nyndg,nzndg)     ! Analysis increment for pt
  REAL :: qvincr(nxndg,nyndg,nzndg)     ! Analysis increment for qv
  REAL :: qcincr(nxndg,nyndg,nzndg)     ! Analysis increment for qc
  REAL :: qrincr(nxndg,nyndg,nzndg)     ! Analysis increment for qr
  REAL :: qiincr(nxndg,nyndg,nzndg)     ! Analysis increment for qi
  REAL :: qsincr(nxndg,nyndg,nzndg)     ! Analysis increment for qs
  REAL :: qhincr(nxndg,nyndg,nzndg)     ! Analysis increment for qh
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL :: timscl,dmid
  REAL :: timsum = 0.
  INTEGER :: icall = 0
  SAVE timsum,icall
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Compute scale factor for this application.
!  nudgopt=1    Constant time function
!  nudgopt=2    Triangular time function, max at mid, zero at ends.
!
!-----------------------------------------------------------------------
!
  icall = icall + 1
  IF(nudgopt == 1) THEN
    timscl=ndscale
  ELSE IF(nudgopt == 2) THEN
    dmid=ABS( (curtim-ndstart) - 0.5*(ndstop-ndstart) )
    timscl=ndscale*AMAX1(0.,(2.-(4.*dmid/(ndstop-ndstart))))
  ELSE
    timscl=0.
  END IF
  timsum=timsum+0.5*timscl
  IF (myproc == 0) WRITE(6,'(a,f9.4,a,f9.4,a,f9.4)')                    &
        ' Timeweight:',timscl,'  Accum wgt:',timsum,                    &
        ' Accum Target:',ndgain

  IF (nudgu == 1 .OR. nudgu == 2 ) THEN
    CALL nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx,1,ny-1,2,nz-1,u,uincr,timscl)
  ELSE IF (nudgu == 3 .AND. icall <= 2 ) THEN
    CALL nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx,1,ny-1,2,nz-1,u,uincr,1.0)
  END IF

  IF (nudgv == 1 .OR. nudgv == 2 ) THEN
    CALL nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx-1,1,ny,2,nz-1,v,vincr,timscl)
  ELSE IF (nudgv == 3 .AND. icall <= 2 ) THEN
    CALL nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx-1,1,ny,2,nz-1,v,vincr,1.0)
  END IF

  IF (nudgw == 1 .OR. nudgw == 2 ) THEN
    CALL nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx-1,1,ny-1,2,nz-1,w,wincr,timscl)
  ELSE IF (nudgw == 3 .AND. icall <= 2 ) THEN
    CALL nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx-1,1,ny-1,2,nz-1,w,wincr,1.0)
  END IF

  IF (nudgp == 1 .OR. nudgp == 2 ) THEN
    CALL nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx-1,1,ny-1,2,nz-1,pprt,pincr,timscl)
  ELSE IF (nudgp == 3 .AND. icall <= 2 ) THEN
    CALL nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx-1,1,ny-1,2,nz-1,pprt,pincr,1.0)
  END IF

  IF (nudgpt == 1 .OR. nudgpt == 2 ) THEN
    CALL nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx-1,1,ny-1,2,nz-1,ptprt,ptincr,timscl)
  ELSE IF (nudgpt == 3 .AND. icall <= 2 ) THEN
    CALL nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx-1,1,ny-1,2,nz-1,ptprt,ptincr,1.0)
  END IF

  IF (nudgqv == 1 .OR. nudgqv == 2 ) THEN
    CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx-1,1,ny-1,2,nz-1,qv,qvincr,timscl)
  ELSE IF (nudgqv == 3 .AND. icall <= 2 ) THEN
    CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  1,nx-1,1,ny-1,2,nz-1,qv,qvincr,1.0)
  END IF

  IF (P_QC > 0) THEN
    IF (nudgqc == 1 .OR. nudgqc == 2 ) THEN
      CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                         &
                1,nx-1,1,ny-1,2,nz-1,qscalar(:,:,:,P_QC),qcincr,timscl)
    ELSE IF (nudgqc == 3 .AND. icall <= 2 ) THEN
      CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                         &
                1,nx-1,1,ny-1,2,nz-1,qscalar(:,:,:,P_QC),qcincr,1.0)
    END IF
  END IF

  IF (P_QR > 0) THEN
    IF (nudgqr == 1 .OR. nudgqr == 2 ) THEN
      CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                    1,nx-1,1,ny-1,2,nz-1,qscalar(:,:,:,P_QR),qrincr,timscl)
    ELSE IF (nudgqr == 3 .AND. icall <= 2 ) THEN
      CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                    1,nx-1,1,ny-1,2,nz-1,qscalar(:,:,:,P_QR),qrincr,1.0)
    END IF
  END IF
  
  IF (P_QI > 0) THEN
    IF (nudgqi == 1 .OR. nudgqi == 2 ) THEN
      CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                    1,nx-1,1,ny-1,2,nz-1,qscalar(:,:,:,P_QI),qiincr,timscl)
    ELSE IF (nudgqi == 3 .AND. icall <= 2 ) THEN
      CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                    1,nx-1,1,ny-1,2,nz-1,qscalar(:,:,:,P_QI),qiincr,1.0)
    END IF
  END IF
  
  IF (P_QS > 0) THEN
    IF (nudgqs == 1 .OR. nudgqs == 2 ) THEN
      CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                    1,nx-1,1,ny-1,2,nz-1,qscalar(:,:,:,P_QS),qsincr,timscl)
    ELSE IF (nudgqs == 3 .AND. icall <= 2 ) THEN
      CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                    1,nx-1,1,ny-1,2,nz-1,qscalar(:,:,:,P_QS),qsincr,1.0)
    END IF
  END IF
  
  IF (P_QH > 0) THEN
    IF (nudgqh == 1 .OR. nudgqh == 2 ) THEN
      CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                    1,nx-1,1,ny-1,2,nz-1,qscalar(:,:,:,P_QH),qhincr,timscl)
    ELSE IF (nudgqh == 3 .AND. icall <= 2 ) THEN
      CALL nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                    1,nx-1,1,ny-1,2,nz-1,qscalar(:,:,:,P_QH),qhincr,1.0)
    END IF
  END IF

  RETURN
END SUBROUTINE nudgeall
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE NUDGEVAR                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE nudgevar(nx,ny,nz,nxndg,nyndg,nzndg,                         &
           ibeg,iend,jbeg,jend,kbeg,kend,                               &
           var,varincr,timscl)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  March, 1998
!
!  MODIFICATION HISTORY:
!
!  07/10/2001 (K. Brewster)
!  Added increment array dimensions to argument list for consistency
!  in array dimension statements.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nxndg    Number x grid points for IAU
!    nyndg    Number y grid points for IAU
!    nzndg    Number z grid points for IAU
!
!    var      Variable to be nudged
!    varincr  Increment to apply to variable over time
!    timscl   Scale factor to determine increment to apply
!             at this time
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: nxndg,nyndg,nzndg
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
  REAL :: var(nx,ny,nz)
  REAL :: varincr(nxndg,nyndg,nzndg)
  REAL :: timscl
!
!-----------------------------------------------------------------------
!
!  Misc local variables
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
  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        var(i,j,k)=var(i,j,k)+timscl*varincr(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE nudgevar
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE NUDGEPOS                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE nudgepos(nx,ny,nz,nxndg,nyndg,nzndg,                         &
           ibeg,iend,jbeg,jend,kbeg,kend,                               &
           var,varincr,timscl)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  March, 1998
!
!  MODIFICATION HISTORY:
!  07/10/2001 (K. Brewster)
!  Added increment array dimensions to argument list for consistency
!  in array dimension statements.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nxndg    Number x grid points for IAU
!    nyndg    Number y grid points for IAU
!    nzndg    Number z grid points for IAU
!
!    var      Variable to be nudged
!    varincr  Increment to apply to variable over time
!    timscl   Scale factor to determine increment to apply
!             at this time
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: nxndg,nyndg,nzndg
  REAL :: var(nx,ny,nz)
  REAL :: varincr(nxndg,nyndg,nzndg)
  REAL :: timscl
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        var(i,j,k)=MAX(0.,(var(i,j,k)+timscl*varincr(i,j,k)))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE nudgepos
