!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE write_wrf_bdy               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_wrf_bdy(ncid,io_form,times_str,wrfversion,ifile,     &
                         nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,         &
                         nxlg_wrf,nylg_wrf,fzone_wrf,tintv_bdyin,     &
                         ubdy3dtemp1,ubdy3dtemp2,                     &
                         vbdy3dtemp1,vbdy3dtemp2,                     &
                         tbdy3dtemp1,tbdy3dtemp2,                     &
                         pbdy3dtemp1,pbdy3dtemp2,                     &
                         qbdy3dtemp1,qbdy3dtemp2,                     &
                         mbdy2dtemp1,mbdy2dtemp2,                     &
                         bdys,bdyn,bdyw,bdye,                         &
                         blgw,blge,blgs,blgn,                         &
                         temx,temy,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!   Write WRF boundary file
!
!-----------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(IN)  :: io_form
  CHARACTER(LEN=*), INTENT(IN) :: times_str
  REAL,    INTENT(IN)  :: wrfversion
  INTEGER, INTENT(IN)  :: ifile
  INTEGER, INTENT(IN)  :: nx_wrf, ny_wrf, nz_wrf
  INTEGER, INTENT(IN)  :: spec_bdy_width
  INTEGER, INTENT(IN)  :: nxlg_wrf, nylg_wrf
  INTEGER, INTENT(IN)  :: fzone_wrf
  INTEGER, INTENT(IN)  :: tintv_bdyin
  REAL,    INTENT(IN)  :: ubdy3dtemp1(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: ubdy3dtemp2(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: vbdy3dtemp1(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: vbdy3dtemp2(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: tbdy3dtemp1(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: tbdy3dtemp2(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: pbdy3dtemp1(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: pbdy3dtemp2(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: qbdy3dtemp1(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: qbdy3dtemp2(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(IN)  :: mbdy2dtemp1(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: mbdy2dtemp2(nx_wrf,ny_wrf)

  REAL,    INTENT(OUT) :: bdys(nx_wrf,  nz_wrf,spec_bdy_width)
  REAL,    INTENT(OUT) :: bdyn(nx_wrf,  nz_wrf,spec_bdy_width)
  REAL,    INTENT(OUT) :: bdyw(ny_wrf,  nz_wrf,spec_bdy_width)
  REAL,    INTENT(OUT) :: bdye(ny_wrf,  nz_wrf,spec_bdy_width)
  REAL,    INTENT(OUT) :: blgs(nxlg_wrf,nz_wrf,spec_bdy_width)
  REAL,    INTENT(OUT) :: blgn(nxlg_wrf,nz_wrf,spec_bdy_width)
  REAL,    INTENT(OUT) :: blgw(nylg_wrf,nz_wrf,spec_bdy_width)
  REAL,    INTENT(OUT) :: blge(nylg_wrf,nz_wrf,spec_bdy_width)

  REAL,    INTENT(OUT) :: temx(nx_wrf,nz_wrf,spec_bdy_width)
  REAL,    INTENT(OUT) :: temy(ny_wrf,nz_wrf,spec_bdy_width)

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  REAL :: bdyintv

  INTEGER :: nq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  bdyintv = REAL(tintv_bdyin)

  ! U wind
  CALL stuff_bdy(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,'U',               &
                ubdy3dtemp1,bdys,bdyn,bdyw,bdye)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'U_B','X',           &
                'bdy x-wind component',bdyw,bdye,bdys,bdyn,             &
                nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,fzone_wrf,          &
                blgw,blge,blgs,blgn,nxlg_wrf,nylg_wrf-1,nz_wrf-1,       &
                temx,temy)

  CALL stuff_bdytend(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,'U',           &
                ubdy3dtemp1,ubdy3dtemp2,bdyintv,bdys,bdyn,bdyw,bdye)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'U_BT','X',          &
                'bdy tend x-wind component',bdyw,bdye,bdys,bdyn,        &
                nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,fzone_wrf,          &
                blgw,blge,blgs,blgn,nxlg_wrf,nylg_wrf-1,nz_wrf-1,       &
                temx,temy)

  ! V wind
  CALL stuff_bdy(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,'V',               &
                vbdy3dtemp1,bdys,bdyn,bdyw,bdye)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'V_B','Y',           &
                'bdy y-wind component',                                 &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf,nz_wrf-1,temx,temy)

  CALL stuff_bdytend(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,'V',           &
                vbdy3dtemp1,vbdy3dtemp2,bdyintv,bdys,bdyn,bdyw,bdye)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'V_BT','Y',          &
                'bdy tend y-wind component',                            &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf,nz_wrf-1,temx,temy)

  bdyw(:,:,:) = 0.0
  bdye(:,:,:) = 0.0
  bdys(:,:,:) = 0.0
  bdyn(:,:,:) = 0.0

  ! Z wind
  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'W_B','Z',           &
                'bdy z-wind component',                                 &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,nz_wrf,temx,temy)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'W_BT','Z',          &
                'bdy tend z-wind component',                            &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,nz_wrf,temx,temy)

  ! Geopotential
  CALL stuff_bdy(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,'H',               &
                pbdy3dtemp1,bdys,bdyn,bdyw,bdye)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'PH_B','Z',          &
                'bdy perturbation geopotential',                        &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,nz_wrf,temx,temy)

  CALL stuff_bdytend(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,'H',           &
                pbdy3dtemp1,pbdy3dtemp2,bdyintv,bdys,bdyn,bdyw,bdye)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'PH_BT','Z',         &
                'bdy tend perturbation geopotential',                   &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,nz_wrf,temx,temy)

  ! Potential temperature perturbation
  CALL stuff_bdy(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,'T',               &
                 tbdy3dtemp1,bdys,bdyn,bdyw,bdye)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'T_B','',            &
                'bdy perturbation potential temperature (theta-t0)',    &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,nz_wrf-1,temx,temy)

  CALL stuff_bdytend(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,'T',           &
                tbdy3dtemp1,tbdy3dtemp2,bdyintv,bdys,bdyn,bdyw,bdye)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'T_BT','',           &
               'bdy tend perturbation potential temperature (theta-t0)',&
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,nz_wrf-1,temx,temy)

  ! mu
  CALL stuff_bdy2d(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,                 &
                mbdy2dtemp1,bdys,bdyn,bdyw,bdye)

  CALL writebdy2d(ncid,io_form,times_str(1:19),ifile,'MU_B','',         &
                'bdy perturbation dry air mass in column',              &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,temx,temy)

  CALL stuff_bdytend2d(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,             &
                mbdy2dtemp1,mbdy2dtemp2,bdyintv,bdys,bdyn,bdyw,bdye)

  CALL writebdy2d(ncid,io_form,times_str(1:19),ifile,'MU_BT','',        &
                'bdy tend perturbation dry air mass in column',         &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,temx,temy)

  ! Water vapor mixing ratio
  CALL stuff_bdy(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,'T',               &
                qbdy3dtemp1,bdys,bdyn,bdyw,bdye)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'QVAPOR_B','',       &
                'Water vapor mixing ratio',                             &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,nz_wrf-1,temx,temy)

  CALL stuff_bdytend(nx_wrf,ny_wrf,nz_wrf,spec_bdy_width,'T',           &
                qbdy3dtemp1,qbdy3dtemp2,bdyintv,bdys,bdyn,bdyw,bdye)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'QVAPOR_BT','',      &
                'Water vapor mixing ratio',                             &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,nz_wrf-1,temx,temy)

  !
  ! Although those variables were dumped in wrfbdy_d01, they
  ! are not set to any meaningful value in real_em.f
  ! those are place holders for nested runs.
  !

  bdyw(:,:,:) = 0.0
  bdye(:,:,:) = 0.0
  bdys(:,:,:) = 0.0
  bdyn(:,:,:) = 0.0

  DO nq = 1,nscalar_wrf
    CALL writebdy(ncid,io_form,times_str(1:19),ifile,                   &
                  TRIM(qnames_wrf(nq))//'_B','',TRIM(qdescp_wrf(nq)),   &
                  bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,             &
                  spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,         &
                  nxlg_wrf-1,nylg_wrf-1,nz_wrf-1,temx,temy)

    CALL writebdy(ncid,io_form,times_str(1:19),ifile,                   &
                  TRIM(qnames_wrf(nq))//'_BT','',TRIM(qdescp_wrf(nq)),  &
                  bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,             &
                  spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,         &
                  nxlg_wrf-1,nylg_wrf-1,nz_wrf-1,temx,temy)
  END DO

  IF (P_QNI_wrf > 0) THEN
  ! QNICE
  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'QNICE_B',           &
                '','Ice Number concentration',                          &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,nz_wrf-1,temx,temy)

  CALL writebdy(ncid,io_form,times_str(1:19),ifile,'QNICE_BT',          &
                '','Ice Number concentration',                          &
                bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,               &
                spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,           &
                nxlg_wrf-1,nylg_wrf-1,nz_wrf-1,temx,temy)
  END IF

  IF (wrfversion >= 3.1) THEN

    CALL writebdy2d(ncid,io_form,times_str(1:19),ifile,'HT_SHAD_B','',  &
                  'bdy Height of orographic shadow',                    &
                  bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,             &
                  spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,         &
                  nxlg_wrf-1,nylg_wrf-1,temx,temy)

    CALL writebdy2d(ncid,io_form,times_str(1:19),ifile,'HT_SHAD_BT','', &
                  'bdy tend Height of orographic shadow',               &
                  bdyw,bdye,bdys,bdyn,nx_wrf,ny_wrf,nz_wrf,             &
                  spec_bdy_width,fzone_wrf,blgw,blge,blgs,blgn,         &
                  nxlg_wrf-1,nylg_wrf-1,temx,temy)
  END IF

  istatus = 0

  RETURN
END SUBROUTINE write_wrf_bdy
