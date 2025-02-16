
SUBROUTINE join_hdf (fileheader,nxsm,nysm,nz,nzsoil,nstyps,nxlg,nylg,  &
                     buf_r,buf_rsm,buf_rsoil,buf_rsmsoil,              &
                     buf_i,buf_ism,buf_i16,buf_i16sm,                  &
                     buf_i16soil,buf_i16soilsm,                        &
                     buf_r1,buf_r2,sstat)

  IMPLICIT NONE

  INCLUDE 'mp.inc'
  INCLUDE 'hdf.f90'    ! HDF4 library include file

  CHARACTER(LEN=*), INTENT(IN) :: fileheader

  INTEGER, INTENT(IN)  :: nxsm,nysm,nz,nzsoil,nstyps,nxlg,nylg
  REAL,    INTENT(OUT) :: buf_r(nxlg,nylg,nz), buf_rsm(nxsm,nysm,nz)
  REAL,    INTENT(OUT) :: buf_rsoil(nxlg,nylg,nzsoil,0:nstyps),      &
                          buf_rsmsoil(nxsm,nysm,nzsoil,0:nstyps)
  INTEGER, INTENT(OUT) :: buf_i(nxlg,nylg,nz), buf_ism(nxsm,nysm,nz)
  INTEGER(KIND=selected_int_kind(4)), INTENT(OUT) ::                 &
                                        buf_i16(nxlg,nylg,nz),       &
                                        buf_i16sm(nxsm,nysm,nz)
  INTEGER(KIND=selected_int_kind(4)), INTENT(OUT) ::                 &
                       buf_i16soil(nxlg,nylg,nzsoil,0:nstyps),       &
                       buf_i16soilsm(nxsm,nysm,nzsoil,0:nstyps)
  REAL,    INTENT(OUT) :: buf_r1(nxsm+nysm+nz), buf_r2(nxlg+nylg+nz)
                          ! 1D arrays or max and min
  INTEGER, INTENT(OUT) :: sstat
                          ! join status: sstat=1 if error encountered

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=256) :: filename
  INTEGER :: fi, fj, i, j, k, l
  INTEGER :: kk
  INTEGER :: nxin, nyin, nzin, nzsoilin

  REAL    :: amax, amin, scale

  INTEGER :: ierr

  INTEGER :: sd_id,sd_id1,sd_id2,ndata,nattr,istat,aindex,dindex
  INTEGER :: sds_id,sds_id1,sds_id2,rank,dims(6),dtype,ndattr,adtype
  INTEGER :: sds_index

  CHARACTER (LEN=256 ) :: name, aname
  CHARACTER (LEN=1024) :: char_attr
  CHARACTER (LEN=1), ALLOCATABLE :: lchar_attr(:)
  INTEGER :: nvalues

  INTEGER :: size(4),start(4),stride(4),strideout(4),sizeout(4)
  INTEGER :: x_off,y_off,i0,j0

  INTEGER :: comp_code, comp_prm(1)

  INTEGER :: itmp

! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Version 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver410a,fmtver500a
  CHARACTER (LEN=40) :: fmtver410b,fmtver500b
  INTEGER  :: intver,intver410,intver500

  PARAMETER (fmtver410a='* 004.10 GrADS Soilvar Data',intver410=410)
  PARAMETER (fmtver500a='* 005.00 GrADS Soilvar Data',intver500=500)
  PARAMETER (fmtver410b='004.10 HDF4 Coded Data')
  PARAMETER (fmtver500b='005.00 HDF4 Coded Data')

  CHARACTER (LEN=40) :: fmtverin

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  INTEGER :: sfcreate, sfrdata, sfrnatt, sfscompress, sfselect,        &
             sfsnatt, sfwdata, sffattr, sfendacc, sfn2index

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  nxlg = (nxsm-3)*nproc_x+3
!  nylg = (nysm-3)*nproc_y+3

  stride = 1
  start = 0

  sstat = 0

!-----------------------------------------------------------------------
!
!  Open file 0101, read in its attributes and check the
!  dimensions in the file against the dimensions passed in.
!  Open the joined file as well.
!
!-----------------------------------------------------------------------

  CALL hdfopen(trim(fileheader)//'_001001',1,sd_id1)
  IF (sd_id1 < 0) THEN
    WRITE (6,*) "JOIN_HDF: ERROR opening ",                              &
                 trim(fileheader)//'_0101'," for reading."
    sstat = 1
    RETURN
  END IF

  CALL hdfopen(trim(fileheader),2,sd_id2)
  IF (sd_id2 < 0) THEN
    WRITE (6,*) "JOIN_HDF: ERROR creating HDF4 file: ",                  &
                trim(fileheader)
    sstat = 1
    GO TO 600
  END IF

  CALL hdfinfo(sd_id1,ndata,nattr,istat)

  fmtverin = ""
  CALL hdfrdc(sd_id1,40,"fmtver",fmtverin,istat)

! The following code is a dangerous practice
! but it may be the only way to distinguish
! versions prior to 500.
!
!  IF (fmtverin == fmtver500a .OR. fmtverin == fmtver500b) THEN
    intver=intver500
!  ELSE
!    intver=intver410  ! prior to 500, there is no fmtver variable
!    istat=0
!    WRITE(6,'(/1x,a/,a/)')  'WARNING: Incoming data format are ' //    &
!            'older than version 5.00!!! ',                             &
!            'It is to be read as if it version 4.10!!! '
!  END IF

  CALL hdfrdi(sd_id1,"nx",nxin,istat)
  CALL hdfrdi(sd_id1,"ny",nyin,istat)
  IF (nz > 1) THEN
    CALL hdfrdi(sd_id1,"nz",nzin,istat)
  ELSE
    nzin = nz  ! the file is 2-D so it won't have nz
  END IF

  IF (nzsoil > 1) THEN ! Soil levels expected.
    IF (intver >= intver500) THEN ! New version with OUSoil model.
      CALL hdfrdi(sd_id1,"nzsoil",nzsoilin,istat)
    ELSE
      nzsoilin = 2  ! prior to version 500, it is 2 layer soil
    END IF
  ELSE
    nzsoilin = nzsoil
  END IF

  IF ((nxin /= nxsm).OR.(nyin /= nysm).OR.(nzin /= nz).OR.            &
      (nzsoilin /= nzsoil)) THEN
    WRITE (*,*) "ERROR:  mismatch in sizes."
    WRITE (*,*) "nxin,nyin,nzin,nzsoilin: ",nxin,nyin,nzin,nzsoilin
    WRITE (*,*) "nxsm,nysm,nz,nzsoil:   ",nxsm,nysm,nz,nzsoil
    sstat = 1
    RETURN
  END IF

!-----------------------------------------------------------------------
!
!  Read/write header info.
!
!-----------------------------------------------------------------------

  DO aindex = 0,nattr-1
    CALL hdfainfo(sd_id1,aindex,name,dtype,nvalues,istat)
    IF (dtype == dfnt_char8) THEN
      IF (nvalues > 1024) THEN
        ALLOCATE (lchar_attr(nvalues))
        CALL hdfrdc(sd_id1,nvalues,name,lchar_attr,istat)
        CALL hdfwrtc(sd_id2,nvalues,name,lchar_attr,istat)
        DEALLOCATE(lchar_attr)
      ELSE
        CALL hdfrdc(sd_id1,nvalues,name,char_attr,istat)
        CALL hdfwrtc(sd_id2,nvalues,name,char_attr,istat)
      END IF
    ELSE IF (dtype == dfnt_float32) THEN
      istat = sfrnatt(sd_id1, aindex, buf_r)
      IF (istat /= 0) THEN
        WRITE (6,*) "JOIN_HDF: ERROR reading attribute ",trim(name)
        sstat = 1
        GOTO 600
      END IF
      istat = sfsnatt(sd_id2, trim(name), dfnt_float32, nvalues, buf_r)
      IF (istat /= 0) THEN
        WRITE (6,*) "JOIN_HDF: ERROR, writing attribute ",trim(name)
        sstat = 1
        GOTO 600
      END IF
    ELSE IF (dtype == dfnt_int32) THEN
      IF (trim(name) == 'nx') THEN
        CALL hdfwrti(sd_id2, 'nx', nxlg, istat)
      ELSE IF (trim(name) == 'ny') THEN
        CALL hdfwrti(sd_id2, 'ny', nylg, istat)
      ELSE
        istat = sfrnatt(sd_id1, aindex, buf_i)
        IF (istat /= 0) THEN
          WRITE (6,*) "JOIN_HDF: ERROR reading attribute ",trim(name)
          sstat = 1
          GOTO 600
        END IF
        istat = sfsnatt(sd_id2, trim(name), dfnt_int32, nvalues, buf_i)
        IF (istat /= 0) THEN
          WRITE (6,*) "JOIN_HDF: ERROR, writing attribute ",trim(name)
          sstat = 1
          GOTO 600
        END IF
      END IF
    ELSE
      WRITE (6,*) "JOIN_HDF: ERROR, unknown data type for ",          &
         "attribute ", trim(name)
      sstat = 1
      GOTO 600
    END IF
  END DO

!-----------------------------------------------------------------------
!
!  Read/write each data set.
!
!-----------------------------------------------------------------------

  DO dindex = 0,ndata-1
    sds_id1 = sfselect(sd_id1,dindex)

    ! data set

    CALL hdfdinfo(sds_id1,name,rank,dims,dtype,ndattr,istat)

    IF (rank /= 1) THEN   ! x, y, z etc. 1d arrays do not have comp_prm
                          ! and hdf_comp_code  - WYH - why?
      CALL hdfrdi(sds_id1,"hdf_comp_prm",comp_prm,istat)
      IF (comp_prm(1) /= 0)                                             &
        CALL hdfrdi(sds_id1,"hdf_comp_code",comp_code,istat)
    END IF

    DO fj=1,nproc_y
      DO fi=1,nproc_x

        CALL gtsplitfn(fileheader,nproc_x,nproc_y,1,1,fi,fj,            &
                       0,0,1,2,filename,istat)

        CALL hdfopen(filename,1,sd_id)
        IF (sd_id < 0) THEN
          WRITE (6,*) "JOIN_HDF: ERROR opening HDF4 file: ",            &
                      trim(filename)
          sstat = 1
          GO TO 600
        END IF

        sds_index = sfn2index(sd_id, trim(name))
        IF (sds_index == -1) THEN
          WRITE (6,*) "JOIN_HDF: ERROR, variable ",                     &
                       trim(name)," not found in file",trim(filename),"."
          sstat = 1
          GO TO 600
        END IF
        sds_id = sfselect(sd_id, sds_index)

        x_off = (fi-1) * (nxsm-3)
        y_off = (fj-1) * (nysm-3)

        i0 = min(2,fi)
        j0 = min(2,fj)

        size(1) = nxsm
        size(2) = nysm
        size(3) = dims(3)
        size(4) = dims(4) ! EMK Test

        sizeout(1) = nxlg
        sizeout(2) = nylg
        sizeout(3) = dims(3)
!        sizeout(4) = dims(4) ! EMK Test
        sizeout(4) = nstyps+1

        IF (rank == 1) THEN

          IF (dtype /= dfnt_float32) THEN
            WRITE (6,*) "JOIN_HDF: ERROR, unsupported data type  ",    &
                        "for 1-d variable ",trim(name)
            sstat = 1
            GOTO 600
          END IF

          IF (trim(name) == 'x') THEN       ! x
            IF (fj == 1) THEN
              size(1) = nxsm
              istat = sfrdata(sds_id, start, stride, size, buf_r1)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading ",trim(name),"."
                sstat = 1
                GOTO 600
              END IF
              DO i=i0,nxsm
                buf_r2(i+x_off) = buf_r1(i)
              END DO
            END IF
          ELSE IF (trim(name) == 'y') THEN  ! y
            IF (fi == 1) THEN
              size(1) = nysm
              istat = sfrdata(sds_id, start, stride, size, buf_r1)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading ",trim(name),"."
                sstat = 1
                GOTO 600
              END IF
              DO j=j0,nysm
                buf_r2(j+y_off) = buf_r1(j)
              END DO
            END IF
          ELSE IF (trim(name) == 'z') THEN  ! z
            IF (fi == 1 .and. fj == 1) THEN
              size(1) = nz
              istat = sfrdata(sds_id, start, stride, dims, buf_r1)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading ",trim(name),"."
                sstat = 1
                GOTO 600
              END IF
            END IF
          ELSE
            IF (fi == 1 .and. fj == 1) THEN
              istat = sfrdata(sds_id, start, stride, dims, buf_r)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading ",trim(name),"."
                sstat = 1
                GOTO 600
              END IF
            END IF
          END IF

        ELSE

          IF (dtype == dfnt_float32) THEN

            IF(intver <= intver410                .AND.                 &
               (TRIM(name) == "tsfc"   .OR. TRIM(name) == "tsoil"  .OR. &
                TRIM(name) == "wetsfc" .OR. TRIM(name) == "wetdp")) THEN
              istat = sfrdata(sds_id, start, stride, size, buf_rsm)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),"."
                sstat = 1
                GOTO 600
              END IF
              DO j=j0,nysm
                DO i=i0,nxsm
                  buf_r(i+x_off,j+y_off,1) = buf_rsm(i,j,1)
                END DO
              END DO
            ELSE IF(intver >= intver500           .AND.                &
                    (TRIM(name) == "tsoil" .OR. TRIM(name) == "qsoil")) THEN
              istat = sfrdata(sds_id, start, stride, size, buf_rsmsoil)
              IF (istat /= 0) THEN
                WRITE (6,*) "SPLIT_HDF: ERROR reading ",trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF
              DO kk=0,nstyps
                DO k=1,nzsoil
                  DO j=j0,nysm
                    DO i=i0,nxsm
                      buf_rsoil(i+x_off,j+y_off,k,kk) =                &
                                                buf_rsmsoil(i,j,k,kk)
                    END DO
                  END DO
                END DO
              END DO
            ELSE IF (intver >= intver500 .AND.                         &
                      (TRIM(name) == "zpsoil")) THEN
              istat = sfrdata(sds_id, start, stride,size,buf_rsmsoil(1,1,1,0))
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading ",trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF
              DO k=1,nzsoil
                DO j=j0,nysm
                  DO i=i0,nxsm
                    buf_rsoil(i+x_off,j+y_off,k,0) = &
                        buf_rsmsoil(i,j,k,0)
                  END DO
                END DO
              END DO
            ELSE ! Not a soil variable

              istat = sfrdata(sds_id, start, stride, size, buf_rsm)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading ",trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF
              DO k=1,nz
                DO j=j0,nysm
                  DO i=i0,nxsm
                    buf_r(i+x_off,j+y_off,k) = buf_rsm(i,j,k)
                  END DO
                END DO
              END DO
            END IF

          ELSE IF (dtype == dfnt_int32) THEN

            istat = sfrdata(sds_id, start, stride, size, buf_ism)
            IF (istat /= 0) THEN
              WRITE (6,*) "JOIN_HDF: ERROR reading ",trim(name),", exiting"
              sstat = 1
              GOTO 600
            END IF
            DO k=1,nz
              DO j=j0,nysm
                DO i=i0,nxsm
                  buf_i(i+x_off,j+y_off,k) = buf_ism(i,j,k)
                END DO
              END DO
            END DO

          ELSE IF (dtype == dfnt_int16) THEN

            IF (rank == 2) THEN

              aindex = sffattr(sds_id, "max")
              istat = sfrnatt(sds_id, aindex, amax)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading max for ",        &
                   trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF
              aindex = sffattr(sds_id, "min")
              istat = sfrnatt(sds_id, aindex, amin)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading min for ",        &
                   trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF
              istat = sfrdata(sds_id, start, stride, size, buf_i16sm)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading ",trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF

              scale = (amax - amin) / 65534.0
              DO j=j0,nysm
                DO i=i0,nxsm
                  buf_r(i+x_off,j+y_off,1) =                           &
                            scale * (buf_i16sm(i,j,1) + 32767) + amin
                END DO
              END DO

            ELSE IF (rank == 3) THEN

              aindex = sffattr(sds_id, "max")
              istat = sfrnatt(sds_id, aindex, buf_r2)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading max for ",        &
                   trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF
              aindex = sffattr(sds_id, "min")
              istat = sfrnatt(sds_id, aindex, buf_r1)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading min for ",        &
                   trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF

              IF( TRIM(name) == "zpsoil") THEN
                istat = sfrdata(sds_id, start, stride, size, buf_i16soilsm(:,:,:,0))
              ELSE
                istat = sfrdata(sds_id, start, stride, size, buf_i16sm)
              END IF
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading ",trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF

              IF( TRIM(name) == "zpsoil" ) THEN
                DO k = 1, nzsoil
                  scale = (buf_r2(k) - buf_r1(k)) / 65534.0
                  DO j=j0,nysm
                    DO i=i0,nxsm
                      buf_rsoil(i+x_off,j+y_off,k,0) = scale *        &
                          (buf_i16soilsm(i,j,k,0) + 32767) + buf_r1(k)
                    END DO
                  END DO
                END DO

              ELSE IF( TRIM(name) == "stypfrct" ) THEN
                DO k=1, nstyps
                  scale = (buf_r2(k) - buf_r1(k)) / 65534.0
                  DO j=j0,nysm
                    DO i=i0,nxsm
                      buf_r(i+x_off,j+y_off,k) =                        &
                           scale * (buf_i16sm(i,j,k) + 32767) + buf_r1(k)
                    END DO
                  END DO
                END DO

              ELSE IF( TRIM(name) == "wetcanp" ) THEN
                DO k=1, nstyps+1
                  scale = (buf_r2(k) - buf_r1(k)) / 65534.0
                  DO j=j0,nysm
                    DO i=i0,nxsm
                      buf_r(i+x_off,j+y_off,k) =                        &
                           scale * (buf_i16sm(i,j,k) + 32767) + buf_r1(k)
                    END DO
                  END DO
                END DO

              ELSE
                DO k=1,nz
                  scale = (buf_r2(k) - buf_r1(k)) / 65534.0
                  DO j=j0,nysm
                    DO i=i0,nxsm
                      buf_r(i+x_off,j+y_off,k) =                        &
                           scale * (buf_i16sm(i,j,k) + 32767) + buf_r1(k)
                    END DO
                  END DO
                END DO
              END IF

            ELSE IF(rank == 4) THEN

              aindex = sffattr(sds_id, "max")
              istat = sfrnatt(sds_id, aindex, buf_r2)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading max for ",        &
                   trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF
              aindex = sffattr(sds_id, "min")
              istat = sfrnatt(sds_id, aindex, buf_r1)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading min for ",        &
                   trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF
              istat = sfrdata(sds_id, start, stride, size, buf_i16soilsm)
              IF (istat /= 0) THEN
                WRITE (6,*) "JOIN_HDF: ERROR reading ",trim(name),", exiting"
                sstat = 1
                GOTO 600
              END IF

              DO l = 0, nstyps
              DO k = 1, nzsoil
                scale = (buf_r2(k) - buf_r1(k)) / 65534.0
                DO j=j0,nysm
                  DO i=i0,nxsm
                    buf_rsoil(i+x_off,j+y_off,k,l) = scale *        &
                        (buf_i16soilsm(i,j,k,l) + 32767) + buf_r1(k)
                  END DO
                END DO
              END DO
              END DO

            ELSE

              WRITE (6,*) "JOIN_HDF: ERROR, unsupported rank,",rank,  &
                 " for 16 bit remapped data"
              GOTO 600

            END IF

          ELSE

            WRITE (6,*) "JOIN_HDF: ERROR, unknown data type for ",    &
               "attribute ", trim(aname), " of sds ",trim(name)
            sstat = 1
            GOTO 600

          END IF

        END IF

        CALL hdfclose(sd_id,istat)
        IF (istat /= 0) THEN
          WRITE (6,*) "JOIN_HDF: ERROR on close of file ",trim(filename)
          sstat = 1
          GOTO 600
        END IF

      END DO ! fi
    END DO ! fj

    IF (rank == 1) THEN

      IF (trim(name) == 'x') THEN       ! x
        sizeout(1) = nxlg
        sds_id2 = sfcreate(sd_id2, trim(name), dfnt_float32, 1, sizeout)
        istat = sfwdata(sds_id2, start, stride, sizeout, buf_r2)
        IF (istat /= 0) THEN
          WRITE (6,*) "JOIN_HDF: ERROR writing ",trim(name),           &
             " to file ",trim(fileheader)," , exiting"
          sstat = 1
          GOTO 600
        END IF
      ELSE IF (trim(name) == 'y') THEN  ! y
        sizeout(1) = nylg
        sds_id2 = sfcreate(sd_id2,trim(name),dfnt_float32,1,sizeout)
        istat = sfwdata(sds_id2, start, stride, sizeout, buf_r2)
        IF (istat /= 0) THEN
          WRITE (6,*) "JOIN_HDF: ERROR writing ",trim(name),           &
             " to file ",trim(fileheader)," , exiting"
          sstat = 1
          GOTO 600
        END IF
      ELSE IF (trim(name) == 'z') THEN  ! z
        sizeout(1) = nz
        sds_id2 = sfcreate(sd_id2,trim(name),dfnt_float32,1,sizeout)
        istat   = sfwdata(sds_id2, start, stride, sizeout, buf_r1)
                  ! note, buf_r1 not buf_r2
        IF (istat /= 0) THEN
          WRITE (6,*) "JOIN_HDF: ERROR writing ",trim(name),           &
             " to file ",trim(fileheader)," , exiting"
          sstat = 1
          GOTO 600
        END IF
      ELSE
        sds_id2 = sfcreate(sd_id2,trim(name),dfnt_float32,1,dims)
        istat = sfwdata(sds_id2, start, stride, dims, buf_r)
        IF (istat /= 0) THEN
          WRITE (6,*) "JOIN_HDF: ERROR writing ",trim(name),           &
             " to file ",trim(fileheader)," , exiting"
          sstat = 1
          GOTO 600
        END IF
      END IF

    ELSE

      IF (dtype == dfnt_float32) THEN

        IF (intver <= intver410         .AND.                         &
            (TRIM(name) == "tsfc"   .OR. TRIM(name) == "tsoil" .OR.   &
             TRIM(name) == "wetsfc" .OR. TRIM(name) == "wetdp")) THEN

          sds_id2 = sfcreate(sd_id2,trim(name),dfnt_float32,rank,sizeout)
          IF (comp_prm(1) > 0) THEN
            istat = sfscompress(sds_id2, comp_code, comp_prm)
          END IF
          istat = sfwdata(sds_id2, start, stride, sizeout, buf_r)
          IF (istat /= 0) THEN
            WRITE (6,*) "JOIN_HDF: ERROR writing ",trim(name),        &
               " to file ",trim(fileheader)," , exiting"
            sstat = 1
            GOTO 600
          END IF
        ELSE IF (intver >= intver500          .AND.                    &
                 (TRIM(name) == "tsoil" .OR. TRIM(name) == "qsoil")) THEN
          sds_id2 = sfcreate(sd_id2,trim(name),dfnt_float32,rank,sizeout)
          IF (comp_prm(1) > 0) THEN
            istat = sfscompress(sds_id2, comp_code, comp_prm)
          END IF
          istat = sfwdata(sds_id2, start, stride, sizeout, buf_rsoil)
          IF (istat /= 0) THEN
            WRITE (6,*) "JOIN_HDF: ERROR writing ",trim(name),         &
               " to file ",trim(fileheader)," , exiting"
            sstat = 1
            GOTO 600
          END IF

        ELSE IF (intver >= intver500 .AND. (TRIM(name) == "zpsoil")) THEN
          sds_id2 = sfcreate(sd_id2,trim(name),dfnt_float32,rank,sizeout)
          IF (comp_prm(1) > 0) THEN
            istat = sfscompress(sds_id2, comp_code, comp_prm)
          END IF
          istat = sfwdata(sds_id2, start, stride, sizeout, buf_rsoil(1,1,1,0))
          IF (istat /= 0) THEN
            WRITE (6,*) "JOIN_HDF: ERROR writing ",trim(name),         &
               " to file ",trim(fileheader)," , exiting"
            sstat = 1
            GOTO 600
          END IF

        ELSE

          sds_id2 = sfcreate(sd_id2,trim(name),dfnt_float32,rank,sizeout)
          IF (comp_prm(1) > 0) THEN
            istat = sfscompress(sds_id2, comp_code, comp_prm)
          END IF
          istat = sfwdata(sds_id2, start, stride, sizeout, buf_r)
          IF (istat /= 0) THEN
            WRITE (6,*) "JOIN_HDF: ERROR writing ",trim(name),           &
               " to file ",trim(fileheader)," , exiting"
            sstat = 1
            GOTO 600
          END IF
        END IF

      ELSE IF (dtype == dfnt_int32) THEN

        sds_id2 = sfcreate(sd_id2,trim(name),dfnt_int32,rank,sizeout)
        IF (comp_prm(1) > 0) THEN
          istat = sfscompress(sds_id2, comp_code, comp_prm)
        END IF
        istat = sfwdata(sds_id2, start, stride, sizeout, buf_i)
        IF (istat /= 0) THEN
          WRITE (6,*) "JOIN_HDF: ERROR writing ",trim(name),          &
             " to file ",trim(fileheader)," , exiting"
          sstat = 1
          GOTO 600
        END IF

      ELSE IF (dtype == dfnt_int16) THEN

        sds_id2 = sfcreate(sd_id2,trim(name),dfnt_int16,rank,sizeout)
        IF (comp_prm(1) > 0) THEN
          istat = sfscompress(sds_id2, comp_code, comp_prm)
        END IF

        IF (rank == 2) THEN

          CALL a3dmax0lcl(buf_r,1,nxlg,1,nxlg,1,nylg,1,nylg,1,1,1,1,amax,amin)
          IF (ABS(amax-amin) < 1.0E-10) THEN
            scale = 65534.0
          ELSE
            scale = 65534.0 / (amax - amin)
          END IF
          DO j=1,nylg
            DO i=1,nxlg
              itmp = nint(scale * (buf_r(i,j,1) - amin)) - 32767
              buf_i16(i,j,1) = itmp
            END DO
          END DO
          buf_r1(1) = amin
          buf_r2(1) = amax
          kk = 1

        ELSE IF (rank == 3) THEN

          IF (TRIM(name) == "zpsoil" ) THEN
            DO k = 1, nzsoil
              CALL a3dmax0lcl(buf_rsoil(1,1,k,0),1,nxlg,1,nxlg,1,nylg,1,nylg,  &
               1,nzsoil,1,1,amax,amin)
              buf_r2(k) = amax
              buf_r1(k) = amin

              IF (ABS(buf_r2(k)-buf_r1(k)) < 1.0E-10) THEN
                scale = 65534.0
              ELSE
                scale = 65534.0 / (buf_r2(k) - buf_r1(k))
              END IF
              DO j=1,nylg
                DO i=1,nxlg
                  itmp = nint(scale * (buf_rsoil(i,j,k,0) - buf_r1(k))) - 32767
                  buf_i16soil(i,j,k,0) = itmp
                END DO
              END DO
            END DO
            kk = nzsoil

          ELSE IF (TRIM(name) == "stypfrct" ) THEN
            DO k = 1, nstyps
              CALL a3dmax0lcl(buf_r(1,1,k),1,nxlg,1,nxlg,1,nylg,1,nylg,  &
               1,nstyps,1,1,amax,amin)
              buf_r2(k) = amax
              buf_r1(k) = amin

              IF (ABS(buf_r2(k)-buf_r1(k)) < 1.0E-10) THEN
                scale = 65534.0
              ELSE
                scale = 65534.0 / (buf_r2(k) - buf_r1(k))
              END IF
              DO j=1,nylg
                DO i=1,nxlg
                  itmp = nint(scale * (buf_r(i,j,k) - buf_r1(k))) - 32767
                  buf_i16(i,j,k) = itmp
                END DO
              END DO
            END DO
            kk = nstyps

          ELSE IF (TRIM(name) == "wetcanp" ) THEN
            DO k = 1, nstyps+1
              CALL a3dmax0lcl(buf_r(1,1,k),1,nxlg,1,nxlg,1,nylg,1,nylg,  &
               1,nstyps+1,1,1,amax,amin)
              buf_r2(k) = amax
              buf_r1(k) = amin

              IF (ABS(buf_r2(k)-buf_r1(k)) < 1.0E-10) THEN
                scale = 65534.0
              ELSE
                scale = 65534.0 / (buf_r2(k) - buf_r1(k))
              END IF
              DO j=1,nylg
                DO i=1,nxlg
                  itmp = nint(scale * (buf_r(i,j,k) - buf_r1(k))) - 32767
                  buf_i16(i,j,k) = itmp
                END DO
              END DO
            END DO
            kk = nstyps + 1

          ELSE
            DO k=1,nz
              CALL a3dmax0lcl(buf_r(1,1,k),1,nxlg,1,nxlg,1,nylg,1,nylg,  &
                 1,nz,1,1,buf_r2(k),buf_r1(k))
              IF (ABS(buf_r2(k)-buf_r1(k)) < 1.0E-10) THEN
                scale = 65534.0
              ELSE
                scale = 65534.0 / (buf_r2(k) - buf_r1(k))
              END IF
              DO j=1,nylg
                DO i=1,nxlg
                  itmp = nint(scale * (buf_r(i,j,k) - buf_r1(k))) - 32767
                  buf_i16(i,j,k) = itmp
                END DO
              END DO
            END DO
            kk = nz
          END IF

        ELSE IF (rank == 4) THEN

          DO l = 0, nstyps
          DO k = 1, nzsoil
            CALL a3dmax0lcl(buf_rsoil(1,1,k,l),1,nxlg,1,nxlg,1,nylg,1,nylg,  &
               1,nzsoil,1,1,amax,amin)
            IF(l == 0) THEN
              buf_r2(k) = amax
              buf_r1(k) = amin
            END IF
            buf_r2(k) = MAX(buf_r2(k),amax)
            buf_r1(k) = MIN(buf_r1(k),amin)
          END DO
          END DO

          DO l = 0, nstyps
          DO k = 1, nzsoil
            IF (ABS(buf_r2(k)-buf_r1(k)) < 1.0E-10) THEN
               scale = 65534.0
            ELSE
               scale = 65534.0 / (buf_r2(k) - buf_r1(k))
            END IF
            DO j=1,nylg
              DO i=1,nxlg
                itmp = nint(scale * (buf_rsoil(i,j,k,l) - buf_r1(k))) - 32767
                buf_i16soil(i,j,k,l) = itmp
              END DO
            END DO
          END DO
          END DO
          kk = nzsoil
        END IF

        IF(rank < 4 .AND. TRIM(name) /= "zpsoil" ) THEN
          istat = sfwdata(sds_id2, start, stride, sizeout, buf_i16)
        ELSE
          istat = sfwdata(sds_id2, start, stride, sizeout, buf_i16soil)
        END IF

        IF (istat /= 0) THEN
          WRITE (6,*) "JOIN_HDF: ERROR writing ",trim(name),  &
             " to file ",trim(fileheader)," , exiting"
          sstat = 1
          GOTO 600
        END IF

      END IF

    END IF

    ! data set attributes

    DO aindex = 0,ndattr-1
      CALL hdfainfo(sds_id1,aindex,aname,adtype,nvalues,istat)

      IF ( dtype == dfnt_int16 .AND. TRIM(aname) == "packed16" ) THEN
        istat = sfsnatt(sds_id2, 'packed16', dfnt_int32, 1, 1)
        CYCLE
      ELSE IF ( dtype == dfnt_int16 .AND. TRIM(aname) == "min" ) THEN
        istat = sfsnatt(sds_id2, 'min', dfnt_float32, kk, buf_r1)
        CYCLE
      ELSE IF ( dtype == dfnt_int16 .AND. TRIM(aname) == "max" ) THEN
        istat = sfsnatt(sds_id2, 'max', dfnt_float32, kk, buf_r2)
        CYCLE
      END IF

      IF (adtype == dfnt_char8) THEN
          IF (nvalues > 1024) THEN
            ALLOCATE (lchar_attr(nvalues))
            CALL hdfrdc(sds_id1,nvalues,aname,lchar_attr,istat)
            CALL hdfwrtc(sds_id2,nvalues,aname,lchar_attr,istat)
            DEALLOCATE(lchar_attr)
          ELSE
            CALL hdfrdc(sds_id1,nvalues,aname,char_attr,istat)
            CALL hdfwrtc(sds_id2,nvalues,aname,char_attr,istat)
          END IF
      ELSE IF (adtype == dfnt_float32) THEN

          IF (size(3) == nzsoil) THEN
            istat = sfrnatt(sds_id1, aindex, buf_rsoil)
            IF (istat /= 0) THEN
              WRITE (6,*) "JOIN_HDF: ERROR reading attribute ",trim(aname)
              sstat = 1
              GOTO 600
            END IF
            istat = sfsnatt(sds_id2, trim(aname), dfnt_float32, nvalues, &
                            buf_rsoil)
            IF (istat /= 0) THEN
              WRITE (6,*) "JOIN_HDF: ERROR, writing attribute ",trim(aname)
              sstat = 1
              GOTO 600
            END IF
          ELSE
            istat = sfrnatt(sds_id1, aindex, buf_r)
            IF (istat /= 0) THEN
              WRITE (6,*) "JOIN_HDF: ERROR reading attribute ",trim(aname)
              sstat = 1
              GOTO 600
            END IF
            istat = sfsnatt(sds_id2, trim(aname), dfnt_float32, nvalues, &
                            buf_r)
            IF (istat /= 0) THEN
              WRITE (6,*) "JOIN_HDF: ERROR, writing attribute ",trim(aname)
              sstat = 1
              GOTO 600
            END IF
          END IF
      ELSE IF (adtype == dfnt_int32) THEN
          istat = sfrnatt(sds_id1, aindex, buf_i)
          IF (istat /= 0) THEN
            WRITE (6,*) "JOIN_HDF: ERROR reading attribute ",trim(aname)
            sstat = 1
            GOTO 600
          END IF
          istat = sfsnatt(sds_id2, trim(aname), dfnt_int32, nvalues, buf_i)
          IF (istat /= 0) THEN
            WRITE (6,*) "JOIN_HDF: ERROR, writing attribute ",trim(aname)
            sstat = 1
            GOTO 600
          END IF
      ELSE
          WRITE (6,*) "JOIN_HDF: ERROR, unknown data type for ", &
             "attribute ", trim(aname)
          sstat = 1
          GOTO 600

      END IF   ! adtype ==

    END DO

    istat = sfendacc(sds_id2)
    IF (istat /= 0) THEN
      WRITE (6,*) "JOIN_HDF: ERROR writing variable ",trim(name)
    END IF

  END DO ! dindex

!-----------------------------------------------------------------------
!
!  Close I/O and deallocate.
!
!----------------------------------------------------------------------

  600   CONTINUE

  CALL hdfclose(sd_id1,istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "JOIN_HDF: ERROR on close of file ",trim(fileheader)//'_0101'
    sstat = 1
  END IF
  CALL hdfclose(sd_id2,istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "JOIN_HDF: ERROR on close of file ",trim(fileheader)
    sstat = 1
  END IF

  RETURN

END SUBROUTINE join_hdf

