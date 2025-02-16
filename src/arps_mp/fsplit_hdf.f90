SUBROUTINE split_hdf(fileheader,nxsm,nysm,nz,nzsoil,nstyps, &
                     buf_r,buf_rsoil4d,buf_i,buf_i16,sstat)

  IMPLICIT NONE

  INCLUDE 'mp.inc'
  INCLUDE 'hdf.f90'    ! HDF4 library include file

  CHARACTER (LEN=*) :: fileheader

  INTEGER ::  nxsm,nysm,nz,nzsoil,nstyps
  INTEGER :: sstat     ! split status: sstat=1 if error encountered

  INTEGER :: nxlg, nylg

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=256) :: filename
  INTEGER :: fi, fj
  INTEGER :: nxin, nyin, nzin,nzsoilin

  REAL :: buf_r(nxsm,nysm,nz)
  REAL :: buf_rsoil4d(nxsm,nysm,nzsoil,0:nstyps)
  INTEGER :: buf_i(nxsm,nysm,nz)
  INTEGER (KIND=selected_int_kind(4)):: buf_i16(nxsm,nysm,nz,0:nstyps)
           ! already assume nz > nzsoil, otherwise there will be trouble

  INTEGER :: sd_id,sd_id2,ndata,nattr,istat,aindex,dindex
  INTEGER :: sds_id,sds_id2,rank,dims(6),dtype,ndattr
  CHARACTER (LEN=256 ) :: name, aname
  CHARACTER (LEN=1024) :: char_attr
  CHARACTER (LEN=1), ALLOCATABLE :: lchar_attr(:)
  INTEGER :: nvalues
  INTEGER :: istart, iend, jstart, jend

  INTEGER :: size(4),start(4),stride(4),startout(4)
  INTEGER :: x_off, y_off

  INTEGER :: comp_code, comp_prm(1)

! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Version 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver410a,fmtver500a,fmtver530a
  CHARACTER (LEN=40) :: fmtver410b,fmtver500b,fmtver530b
  INTEGER  :: intver,intver410,intver500,intver530

  PARAMETER (fmtver410a='* 004.10 GrADS Soilvar Data',intver410=410)
  PARAMETER (fmtver500a='* 005.00 GrADS Soilvar Data',intver500=500)
  PARAMETER (fmtver530a='* 005.30 GrADS Soilvar Data',intver530=530)
  PARAMETER (fmtver410b='004.10 HDF4 Coded Data')
  PARAMETER (fmtver500b='005.00 HDF4 Coded Data')
  PARAMETER (fmtver530b='005.30 HDF4 Coded Data')

  CHARACTER (LEN=40) :: fmtverin

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfcreate, sfrdata, sfrnatt, sfscompress, sfselect,  &
     sfsnatt, sfwdata, sfendacc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = (nxsm-3)*nproc_x+3
  nylg = (nysm-3)*nproc_y+3

  stride = 1
  startout = 0

  sstat = 0

!-----------------------------------------------------------------------
!
!  Open the original file, read in its attributes and check the
!  dimensions in the file against the dimensions passed in.
!
!-----------------------------------------------------------------------

  CALL hdfopen(trim(fileheader),1,sd_id)
  IF (sd_id < 0) THEN
    WRITE (6,*) "SPLIT_HDF: ERROR opening ",                              &
                 trim(fileheader)," for reading."
    sstat = 1
    RETURN
  END IF

  CALL hdfinfo(sd_id,ndata,nattr,istat)

  fmtverin = " "
  CALL hdfrdc(sd_id,40,"fmtver",fmtverin,istat)

! The following code is a dangerous practice
! but it may be the only way to distinguish
! versions prior to 500.
!
  IF (fmtverin == fmtver530a .OR. fmtverin == fmtver530b) THEN
    intver=intver530
  ELSE IF (fmtverin == fmtver500a .OR. fmtverin == fmtver500b) THEN
    intver=intver500
  ELSE
    intver=intver410  ! prior to 500, there is no fmtver variable
    istat=0
    WRITE(6,'(/1x,a/,a/,1x,a/,a/)')   &
          'WARNING: Incoming data format are older than version 5.00!!! ', &
          '          It is to be read as if version 4.10!!! ',             &
          'NOTE: Ignore this WARNING for terain data and surface data,',   &
          '       because both have the same format for version 5.00 and 4.10.'
  END IF

  CALL hdfrdi(sd_id,"nx",nxin,istat)
  CALL hdfrdi(sd_id,"ny",nyin,istat)
  IF (nz > 1) THEN
    CALL hdfrdi(sd_id,"nz",nzin,istat)
  ELSE
    nzin = nz          ! the file is 2-D so it won't have nz
  ENDIF

  IF (nzsoil > 1) THEN ! Soil levels expected.
    IF (intver >= intver500) THEN
                       ! New version with OUSoil model.
      CALL hdfrdi(sd_id,"nzsoil",nzsoilin,istat)
    ELSE
      nzsoilin = 2     ! prior to version 500, it is 2 layer soil
    END IF
  ELSE
    nzsoilin = nzsoil
  END IF

  IF ((nxin /= nxlg).OR.(nyin /= nylg).OR.(nzin /= nz).OR. &
      (nzsoilin /= nzsoil)) THEN
    WRITE (*,*) "ERROR:  mismatch in sizes."
    WRITE (*,*) "nxin,nyin,nzin,nzsoilin: ",nxin,nyin,nzin,nzsoilin
    WRITE (*,*) "nxlg,nylg,nz,nzsoil:   ",nxlg,nylg,nz,nzsoil
    sstat = 1
    RETURN
  END IF

  if ( mp_opt > 0 ) then
      istart = loc_x
      iend   = loc_x
      jstart = loc_y
      jend   = loc_y
  else
      istart = 1
      iend   = nproc_x
      jstart = 1
      jend   = nproc_y
  endif

! DO fj=1,nproc_y
!   DO fi=1,nproc_x

  DO fj=jstart,jend
    DO fi=istart,iend

      x_off = (fi-1) * (nxsm-3)
      y_off = (fj-1) * (nysm-3)

!-----------------------------------------------------------------------
!
!  Create each split file.
!
!-----------------------------------------------------------------------

      CALL gtsplitfn(fileheader,1,1,fi,fj,1,1,0,0,0,2,filename,istat)

      CALL hdfopen(filename,2,sd_id2)
      IF (sd_id2 < 0) THEN
        WRITE (6,*) "SPLIT_HDF: ERROR creating HDF4 file: ",                  &
                    trim(filename)
        sstat = 1
        GO TO 600
      END IF

!-----------------------------------------------------------------------
!
!  Read/write header info.
!
!-----------------------------------------------------------------------

      DO aindex = 0,nattr-1
        CALL hdfainfo(sd_id,aindex,name,dtype,nvalues,istat)
        IF (dtype == dfnt_char8) THEN
          IF (nvalues > 1024) THEN
            ALLOCATE (lchar_attr(nvalues))
            CALL hdfrdc(sd_id,nvalues,name,lchar_attr,istat)
            CALL hdfwrtc(sd_id2,nvalues,name,lchar_attr,istat)
            DEALLOCATE(lchar_attr)
          ELSE
            CALL hdfrdc(sd_id,nvalues,name,char_attr,istat)
            CALL hdfwrtc(sd_id2,nvalues,name,char_attr,istat)
          ENDIF
        ELSE IF (dtype == dfnt_float32) THEN
          !EMK:  Add soil stuff here???
          istat = sfrnatt(sd_id, aindex, buf_r)
          IF (istat /= 0) THEN
            WRITE (6,*) "SPLIT_HDF: ERROR, reading attribute ",trim(name)
            sstat = 1
            GOTO 600
          ENDIF
          istat = sfsnatt(sd_id2, trim(name), dfnt_float32, nvalues, buf_r)
          IF (istat /= 0) THEN
            WRITE (6,*) "SPLIT_HDF: ERROR, writing attribute ",trim(name)
            sstat = 1
            GOTO 600
          ENDIF
        ELSE IF (dtype == dfnt_int32) THEN
          IF (trim(name) == 'nx') THEN
            CALL hdfwrti(sd_id2, 'nx', nxsm, istat)
          ELSE IF (trim(name) == 'ny') THEN
            CALL hdfwrti(sd_id2, 'ny', nysm, istat)
          ELSE
            istat = sfrnatt(sd_id, aindex, buf_i)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR, reading attribute ",trim(name)
              sstat = 1
              GOTO 600
            ENDIF
            istat = sfsnatt(sd_id2, trim(name), dfnt_int32, nvalues, buf_i)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR, writing attribute ",trim(name)
              sstat = 1
              GOTO 600
            ENDIF
          ENDIF
        ELSE
          WRITE (6,*) "SPLIT_HDF: ERROR, unknown data type for ", &
             "attribute ", trim(name)
          sstat = 1
          GOTO 600
        ENDIF
      END DO

!-----------------------------------------------------------------------
!
!  Read/write each data set.
!
!-----------------------------------------------------------------------

      DO dindex = 0,ndata-1
        sds_id = sfselect(sd_id,dindex)

        ! data set

        CALL hdfdinfo(sds_id,name,rank,dims,dtype,ndattr,istat)

        start(1) = x_off
        start(2) = y_off
        start(3) = 0
        start(4) = 0 ! EMK Test
        size(1) = nxsm
        size(2) = nysm
        size(3) = dims(3)
!       size(4) = dims(4)  ! EMK Test
        size(4) = nstyps+1 ! KWT Fix

        IF (rank /= 1) THEN         ! x,y,z etc. 1d array do not write comp_prm.
                                    ! - WYH - why?
          CALL hdfrdi(sds_id,"hdf_comp_prm",comp_prm,istat)
          IF(comp_prm(1) /= 0) CALL hdfrdi(sds_id,"hdf_comp_code",comp_code,istat)

        END IF

        IF (rank == 1) THEN

          IF (dtype /= dfnt_float32) THEN
            WRITE (6,*) "SPLIT_HDF: ERROR, unsupported data type for 1-d ",  &
               " variable ",trim(name)
            sstat = 1
            GOTO 600
          ENDIF

          IF (trim(name) == 'x') THEN       ! x
            size(1) = nxsm
            istat = sfrdata(sds_id, start, stride, size, buf_r)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
              sstat = 1
              GOTO 600
            ENDIF
            sds_id2 = sfcreate(sd_id2, trim(name), dfnt_float32, 1, size)
            istat = sfwdata(sds_id2, startout, stride, size, buf_r)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR writing ",trim(name),  &
                 " to file ",trim(filename)," , exiting"
              sstat = 1
              GOTO 600
            ENDIF
          ELSE IF (trim(name) == 'y') THEN  ! y
            start(1) = y_off
            size(1) = nysm
            istat = sfrdata(sds_id, start, stride, size, buf_r)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
              sstat = 1
              GOTO 600
            ENDIF
            sds_id2 = sfcreate(sd_id2, trim(name), dfnt_float32, 1, size)
            istat = sfwdata(sds_id2, startout, stride, size, buf_r)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR writing ",trim(name),  &
                 " to file ",trim(filename)," , exiting"
              sstat = 1
              GOTO 600
            ENDIF
          ELSE IF (trim(name) == 'z') THEN  ! z
            start = 0
            size(1) = nz
            istat = sfrdata(sds_id, start, stride, dims, buf_r)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
              sstat = 1
              GOTO 600
            ENDIF
            sds_id2 = sfcreate(sd_id2, trim(name), dfnt_float32, 1, size)
            istat = sfwdata(sds_id2, startout, stride, size, buf_r)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR writing ",trim(name),  &
                 " to file ",trim(filename)," , exiting"
              sstat = 1
              GOTO 600
            ENDIF
          ELSE
            start = 0
            startout = 0
            istat = sfrdata(sds_id, start, stride, dims, buf_r)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
              sstat = 1
              GOTO 600
            ENDIF
            sds_id2 = sfcreate(sd_id2, trim(name), dfnt_float32, 1, dims)
            istat = sfwdata(sds_id2, startout, stride, dims, buf_r)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR writing ",trim(name),  &
                 " to file ",trim(filename)," , exiting"
              sstat = 1
              GOTO 600
            ENDIF
          ENDIF

        ELSE

          IF (dtype == dfnt_float32) THEN

            IF (intver <= intver410 .AND. &
                 (TRIM(name) == "tsfc" .OR. &
                  TRIM(name) == "tsoil" .OR. &
                  TRIM(name) == "wetsfc" .OR. &
                  TRIM(name) == "wetdp")) THEN

              istat = sfrdata(sds_id, start, stride, size, buf_r)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
                sstat = 1
                GOTO 600
              ENDIF
              sds_id2 = sfcreate(sd_id2, trim(name), dfnt_float32, rank, size)
              IF (comp_prm(1) > 0) THEN
                istat = sfscompress(sds_id2, comp_code, comp_prm)
              ENDIF
              istat = sfwdata(sds_id2, startout, stride, size, buf_r)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR writing ",trim(name),  &
                   " to file ",trim(filename)," , exiting"
                sstat = 1
                GOTO 600
              ENDIF

            ELSE IF (intver >= intver500 .AND. &
                      (TRIM(name) == "tsoil" .OR. &
                       TRIM(name) == "qsoil")) THEN

              istat = sfrdata(sds_id, start, stride, size, buf_rsoil4d)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
                sstat = 1
                GOTO 600
              ENDIF
              sds_id2 = sfcreate(sd_id2, trim(name), dfnt_float32, rank, size)
              IF (comp_prm(1) > 0) THEN
                istat = sfscompress(sds_id2, comp_code, comp_prm)
              ENDIF
              istat = sfwdata(sds_id2, startout, stride, size, buf_rsoil4d)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR writing ",trim(name),  &
                   " to file ",trim(filename)," , exiting"
                sstat = 1
                GOTO 600
              ENDIF

            ELSE IF (intver >= intver500 .AND. &
                      (TRIM(name) == "zpsoil")) THEN

              istat = sfrdata(sds_id, start, stride, size, buf_rsoil4d(1,1,1,0))
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
                sstat = 1
                GOTO 600
              ENDIF
              sds_id2 = sfcreate(sd_id2, trim(name), dfnt_float32, rank, size)
              IF (comp_prm(1) > 0) THEN
                istat = sfscompress(sds_id2, comp_code, comp_prm)
              ENDIF
              istat = sfwdata(sds_id2, startout, stride, size, buf_rsoil4d(1,1,1,0))
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR writing ",trim(name),  &
                   " to file ",trim(filename)," , exiting"
                sstat = 1
                GOTO 600
              ENDIF

            ELSE ! Not a soil variable

              istat = sfrdata(sds_id, start, stride, size, buf_r)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
                sstat = 1
                GOTO 600
              ENDIF
              sds_id2 = sfcreate(sd_id2, trim(name), dfnt_float32, rank, size)
              IF (comp_prm(1) > 0) THEN
                istat = sfscompress(sds_id2, comp_code, comp_prm)
              ENDIF
              istat = sfwdata(sds_id2, startout, stride, size, buf_r)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR writing ",trim(name),  &
                   " to file ",trim(filename)," , exiting"
                sstat = 1
                GOTO 600
              ENDIF
            ENDIF

          ELSE IF (dtype == dfnt_int32) THEN

            istat = sfrdata(sds_id, start, stride, size, buf_i)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
              sstat = 1
              GOTO 600
            ENDIF
            sds_id2 = sfcreate(sd_id2, trim(name), dfnt_int32, rank, size)
            IF (comp_prm(1) > 0) THEN
              istat = sfscompress(sds_id2, comp_code, comp_prm)
            ENDIF
            istat = sfwdata(sds_id2, startout, stride, size, buf_i)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR writing ",trim(name),  &
                 " to file ",trim(filename)," , exiting"
              sstat = 1
              GOTO 600
            ENDIF

          ELSE IF (dtype == dfnt_int16) THEN

            istat = sfrdata(sds_id, start, stride, size, buf_i16)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
              sstat = 1
              GOTO 600
            ENDIF
            sds_id2 = sfcreate(sd_id2, trim(name), dfnt_int16, rank, size)
            IF (comp_prm(1) > 0) THEN
              istat = sfscompress(sds_id2, comp_code, comp_prm)
            ENDIF
            istat = sfwdata(sds_id2, startout, stride, size, buf_i16)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR writing ",trim(name),  &
                 " to file ",trim(filename)," , exiting"
              sstat = 1
              GOTO 600
            ENDIF

          ELSE

            WRITE (6,*) "SPLIT_HDF: ERROR, unknown data type for ", &
               "attribute ", trim(name)
            sstat = 1
            GOTO 600

          ENDIF

        ENDIF

        ! data set attributes

        DO aindex = 0,ndattr-1
          CALL hdfainfo(sds_id,aindex,aname,dtype,nvalues,istat)
          IF (dtype == dfnt_char8) THEN
            IF (nvalues > 1024) THEN
              ALLOCATE (lchar_attr(nvalues))
              CALL hdfrdc(sds_id,nvalues,aname,lchar_attr,istat)
              CALL hdfwrtc(sds_id2,nvalues,aname,lchar_attr,istat)
              DEALLOCATE(lchar_attr)
            ELSE
              CALL hdfrdc(sds_id,nvalues,aname,char_attr,istat)
              CALL hdfwrtc(sds_id2,nvalues,aname,char_attr,istat)
            ENDIF
          ELSE IF (dtype == dfnt_float32) THEN
            IF (size(3) == nzsoil) THEN
              istat = sfrnatt(sds_id, aindex, buf_rsoil4d)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR, reading attribute ",trim(aname)
                sstat = 1
                GOTO 600
              ENDIF
              istat = sfsnatt(sds_id2, trim(aname), dfnt_float32, nvalues,&
                              buf_rsoil4d)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR, writing attribute ",trim(aname)
                sstat = 1
                GOTO 600
              ENDIF
            ELSE
              istat = sfrnatt(sds_id, aindex, buf_r)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR, reading attribute ",trim(aname)
                sstat = 1
                GOTO 600
              ENDIF
              istat = sfsnatt(sds_id2, trim(aname), dfnt_float32, nvalues, buf_r)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR, writing attribute ",trim(aname)
                sstat = 1
                GOTO 600
              ENDIF
            END IF
          ELSE IF (dtype == dfnt_int32) THEN
            istat = sfrnatt(sds_id, aindex, buf_i)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR, reading attribute ",trim(aname)
              sstat = 1
              GOTO 600
            ENDIF
            istat = sfsnatt(sds_id2, trim(aname), dfnt_int32, nvalues, buf_i)
            IF (istat /= 0) THEN
              WRITE (6,*) "SPLIT_HDF: ERROR, writing attribute ",trim(aname)
              sstat = 1
              GOTO 600
            ENDIF
          ELSE
            WRITE (6,*) "SPLIT_HDF: ERROR, unknown data type for ", &
               "attribute ", trim(aname)
            sstat = 1
            GOTO 600
          ENDIF

        END DO

        istat = sfendacc(sds_id2)
        IF (istat /= 0) THEN
          WRITE (6,*) "SPLIT_HDF: ERROR writing variable ",trim(name)
        END IF

      END DO ! dindex

      CALL hdfclose(sd_id2,istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "SPLIT_HDF: ERROR on close of file ",trim(filename)
        sstat = 1
        GOTO 600
      ENDIF

    END DO ! fi

  END DO ! fj

!-----------------------------------------------------------------------
!
!  Close I/O and return.
!
!----------------------------------------------------------------------

  600   CONTINUE

  CALL hdfclose(sd_id,istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "SPLIT_HDF: ERROR on close of file ",trim(fileheader)
    sstat = 1
  ENDIF

  RETURN

END SUBROUTINE split_hdf

