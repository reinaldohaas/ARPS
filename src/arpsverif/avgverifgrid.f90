!########################################################################
!########################################################################
!######                                                            ######
!######                  PROGRAM AVGVERIFGRID                      ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

PROGRAM avgverifgrid

!-----------------------------------------------------------------------
! 
! PURPOSE:
!
! Reads in HDF files produced by verifgrid, and calculates average
! RMSE and Bias statistics.  Also produces regional statistics.
!
! AUTHOR:  Eric Kemp, May 2000
!
!-----------------------------------------------------------------------
! 
! Use modules
!
!-----------------------------------------------------------------------

  USE verif

!-----------------------------------------------------------------------
!
! Variable declarations
!
!-----------------------------------------------------------------------
  
  IMPLICIT NONE
  
  CHARACTER*132 :: filedir ! temporary

  INTEGER,PARAMETER :: maxntime = 4

  INTEGER :: fdate(6)
  CHARACTER*4 :: amodel,agrid,fmodel,fgrid
  REAL :: adx,fdx
  INTEGER :: fnx,fny,vnx,vny
  INTEGER :: vmapproj
  REAL :: vtrulat1,vtrulat2,vtrulon,vsclfct,vctrlat,vctrlon
  REAL :: corner_lat(2,2),corner_lon(2,2)
  INTEGER :: timesec(maxntime)
  REAL,PARAMETER :: missing = -9999.
  INTEGER :: counter,avgcounter,avgcounterfiles
  REAL :: vdx,vdy
  
  INTEGER :: counter_p(nlevel,maxntime,nvar_p_stats)
  INTEGER :: counter_sfc(maxntime,nvar_sfc_stats)
  REAL,DIMENSION (nlevel,maxntime,nvar_p_stats) :: bias_p,rms_p
  REAL,DIMENSION (maxntime,nvar_sfc_stats) :: bias_sfc,rms_sfc
  REAL,ALLOCATABLE,DIMENSION (:,:,:,:,:) :: diff_p
  REAL,ALLOCATABLE,DIMENSION (:,:,:,:) :: diff_sfc


!..."Master" settings

  CHARACTER*4 :: mamodel,magrid,mfmodel,mfgrid
  REAL :: madx,mfdx
  INTEGER :: mfnx,mfny,mvnx,mvny
  INTEGER :: mvmapproj
  REAL :: mvtrulat1,mvtrulat2,mvtrulon,mvsclfct,mvctrlat,mvctrlon
  REAL :: mcorner_lat(2,2),mcorner_lon(2,2)
  INTEGER :: mtimesec(maxntime)
  INTEGER,ALLOCATABLE :: mcounter(:)
  REAL :: mvdx,mvdy

  INTEGER,ALLOCATABLE :: mcounter_p(:,:,:,:)
  INTEGER,ALLOCATABLE :: mcounter_sfc(:,:,:)
  REAL,ALLOCATABLE,DIMENSION (:,:,:,:) :: mbias_p,mrms_p
  REAL,ALLOCATABLE,DIMENSION (:,:,:) :: mbias_sfc,mrms_sfc
  REAL,ALLOCATABLE,DIMENSION (:,:,:,:,:,:) :: mdiff_p
  REAL,ALLOCATABLE,DIMENSION (:,:,:,:,:) :: mdiff_sfc

  REAL,DIMENSION (nlevel,maxntime,nvar_p_stats) :: avgbias_p,avgrms_p
  INTEGER,DIMENSION(nlevel,maxntime,nvar_p_stats) :: &
                avgbiasfiles_p,avgrmsfiles_p,avgcounterfiles_p
  REAL,DIMENSION (maxntime,nvar_sfc_stats) :: avgbias_sfc,avgrms_sfc
  INTEGER,DIMENSION(maxntime,nvar_sfc_stats) :: &
                avgbiasfiles_sfc,avgrmsfiles_sfc,avgcounterfiles_sfc
  INTEGER :: avgcounter_p(nlevel,maxntime,nvar_p_stats)
  INTEGER :: avgcounter_sfc(maxntime,nvar_sfc_stats)

  INTEGER,ALLOCATABLE :: mfdate(:,:)
!...

  CHARACTER*72 :: line72= &
'------------------------------------------------------------------------'

  CHARACTER*256 :: curfile 
  INTEGER :: i,j,k,l,m

!-----------------------------------------------------------------------
!
! Namelists
! 
!-----------------------------------------------------------------------

  INTEGER :: anx,any
  INTEGER :: numinfiles
  CHARACTER*256 :: infilename(90) ! May need to change this...
  INTEGER :: if_diff = 0
  INTEGER :: numdiffiles
  CHARACTER*256 :: hdiffile(90)

  NAMELIST /input/ anx,any,numinfiles,infilename,if_diff,numdiffiles, &
                   hdiffile

  CHARACTER*256 :: txtout,hdfout
  NAMELIST /output/ txtout,hdfout

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   avgbias_p = 0
   avgbiasfiles_p = 0
   avgrms_p = 0
   avgrmsfiles_p = 0
   avgbias_sfc = 0
   avgbiasfiles_sfc = 0
   avgrms_sfc = 0
   avgrmsfiles_sfc = 0
   avgcounter = 0
   avgcounterfiles = 0
   avgcounter_p = 0
   avgcounterfiles_p = 0
   avgcounter_sfc = 0
   avgcounterfiles_sfc = 0

!-----------------------------------------------------------------------
!
! Read Namelists
! 
!-----------------------------------------------------------------------

  WRITE(6,*) 'Reading NAMELIST input'
  READ(5,input)
  WRITE(6,*) 'Reading NAMELIST output'
  READ(5,output)

  WRITE(6,*) 'anx = ',anx
  WRITE(6,*) 'any = ',any
  WRITE(6,*) 'numinfiles = ',numinfiles

  IF (numinfiles.gt.0) THEN
    DO i = 1,numinfiles
      WRITE(6,*)'i = ',i,' infilename(i) = ',infilename(i) 
    END DO
  END IF

  WRITE(6,*) 'if_diff = ',if_diff
  WRITE(6,*) 'numdiffiles = ',numdiffiles

  IF (numdiffiles.gt.0) THEN
    DO i = 1,numdiffiles
      WRITE(6,*)'i = ',i,' hdiffile(i) = ',hdiffile(i) 
    END DO
  END IF

!-----------------------------------------------------------------------
!
! Allocate arrays and begin processing the data.
! 
!-----------------------------------------------------------------------

  ALLOCATE(diff_p(anx,any,nlevel,maxntime,nvar_sfc_stats), &
           diff_sfc(anx,any,maxntime,nvar_sfc_stats),&
           mdiff_p(anx,any,nlevel,maxntime,nvar_sfc_stats,numinfiles),&
           mdiff_sfc(anx,any,maxntime,nvar_sfc_stats,numinfiles),&
           mcounter_p(nlevel,maxntime,nvar_p_stats,numinfiles),&
           mcounter_sfc(maxntime,nvar_sfc_stats,numinfiles),&
           mbias_p(nlevel,maxntime,nvar_p_stats,numinfiles),&
           mrms_p(nlevel,maxntime,nvar_p_stats,numinfiles),&
           mbias_sfc(maxntime,nvar_sfc_stats,numinfiles),&
           mrms_sfc(maxntime,nvar_sfc_stats,numinfiles), &
           mcounter(numinfiles),mfdate(6,numinfiles))

  DO i = 1,numinfiles

!-----------------------------------------------------------------------
!
!   Read the current file
! 
!-----------------------------------------------------------------------

    curfile = infilename(i)

    CALL rd_verif_stats_diff (curfile, if_diff, hdiffile, fdate, &
                     amodel,agrid,adx,anx,any, &
                     fmodel,fgrid,fdx, &
                     vdx, vdy, vnx, vny, vmapproj, &
                     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat,&
                     vctrlon, corner_lat, corner_lon, &
                     nlevel, maxntime, nvar_p_stats, nvar_sfc_stats, &
                     timesec, pressure, missing, counter, &
                     counter_p, bias_p, rms_p, diff_p, &
                     counter_sfc, bias_sfc, rms_sfc, diff_sfc, &
                     varid_p, varname_p, varunit_p, &
                     varid_sfc, varname_sfc, varunit_sfc)

!-----------------------------------------------------------------------
!
!   If this is the first file, save the various map projection and
!   time information into the "master" arrays.   
! 
!-----------------------------------------------------------------------

    IF (i.EQ.1) THEN
      mamodel = amodel
      magrid = agrid
      mfmodel = fmodel
      mfgrid = fgrid
      madx = adx
      mfdx = fdx
      mfnx = fnx
      mfny = fny
      mvnx = vnx
      mvny = vny
      mvmapproj = vmapproj
      mvtrulat1 = vtrulat1
      mvtrulat2 = vtrulat2
      mvtrulon = vtrulon
      mvsclfct = vsclfct
      mvctrlat = vctrlat
      mvctrlon = vctrlon
      mcorner_lat = corner_lat
      mcorner_lon = corner_lon
      mtimesec = timesec
      mvdx = vdx
      mvdy = vdy      
    END IF
    mfdate(:,i) = fdate

!-----------------------------------------------------------------------
!
!   Compare the current map projection and time information with that
!   stored in the "master" arrays.  If they don't match, abort.
! 
!-----------------------------------------------------------------------

    IF ((mamodel.NE.amodel).OR.(magrid.NE.agrid).OR.&
        (mfmodel.NE.fmodel).OR.(mfgrid.NE.fgrid).OR.&
        (madx.NE.adx).OR.(mfdx.NE.fdx).OR.(mfnx.NE.fnx).OR.&
        (mfnx.NE.fny).OR.(mvnx.NE.vnx).OR.(mvny.NE.vny).OR.&
        (mvmapproj.NE.vmapproj).OR.(mvtrulat1.NE.vtrulat1).OR.&
        (mvtrulat2.NE.vtrulat2).OR.(mvtrulon.NE.vtrulon).OR.&
        (mvsclfct.NE.vsclfct).OR.(mvctrlat.NE.vctrlat).OR.&
        (mvctrlon.NE.mvctrlon).OR.&
        (mcorner_lat(1,1).NE.corner_lat(1,1)).OR.&
        (mcorner_lat(2,1).NE.corner_lat(2,1)).OR.&
        (mcorner_lat(1,2).NE.corner_lat(1,2)).OR.&
        (mcorner_lat(2,2).NE.corner_lat(2,2)).OR.&
        (mcorner_lon(1,1).NE.corner_lon(1,1)).OR.&
        (mcorner_lon(2,1).NE.corner_lon(2,1)).OR.&
        (mcorner_lon(1,2).NE.corner_lon(1,2)).OR.&
        (mcorner_lon(2,2).NE.corner_lon(2,2)).OR.&
        (mvdx.NE.vdx).OR.(mvdy.NE.vdy)) THEN

      WRITE(6,*)'ERROR: Grids do not match!  Aborting...'
      WRITE(6,*)'mamodel = ',TRIM(mamodel),' amodel = ',TRIM(amodel)
      WRITE(6,*)'magrid = ',TRIM(magrid),' agrid = ',TRIM(agrid)
      WRITE(6,*)'mfmodel = ',TRIM(mfmodel),' fmodel = ',TRIM(fmodel)
      WRITE(6,*)'mfgrid = ',TRIM(mfgrid),' fgrid = ',TRIM(fgrid)
      WRITE(6,*)'madx = ',madx,' adx = ',adx
      WRITE(6,*)'mfdx = ',mfdx,' fdx = ',fdx
      WRITE(6,*)'mfnx = ',mfnx,' fnx = ',fnx
      WRITE(6,*)'mfny = ',mfny,' fny = ',fny
      WRITE(6,*)'mvnx = ',mvnx,' vnx = ',vnx
      WRITE(6,*)'mvny = ',mvny,' vny = ',vny
      WRITE(6,*)'mvmapproj = ',mvmapproj,' vmapproj = ', &
                vmapproj
      WRITE(6,*)'mvtrulat1 = ',mvtrulat1,' vtrulat1 = ', &
                vtrulat1
      WRITE(6,*)'mvtrulat2 = ',mvtrulat2,' vtrulat2 = ', &
                vtrulat2
      WRITE(6,*)'mvtrulon = ',mvtrulon,' vtrulon = ',&
                vtrulon
      WRITE(6,*)'mvsclfct = ',mvsclfct,' vsclfct = ',&
                vsclfct
      WRITE(6,*)'mvctrlat = ',mvctrlat,' vctrlat = ',&
                vctrlat
      WRITE(6,*)'mvctrlon = ',mvctrlon,' vctrlon = ',&
                vctrlon
      WRITE(6,*)'mcorner_lat = ',mcorner_lat,' corner_lat = ',&
                corner_lat
      WRITE(6,*)'mcorner_lon = ',mcorner_lon,' corner_lon = ',&
                corner_lon
      WRITE(6,*)'mtimesec = ',mtimesec,' timesec = ',timesec
      WRITE(6,*)'mvdx = ',mvdx,' vdx = ',vdx
      WRITE(6,*)'mvdy = ',mvdy,' vdy = ',vdy
      STOP
    ENDIF
    DO j = 1,maxntime
      IF (mtimesec(j).NE.timesec(j)) THEN
        WRITE(6,*)'ERROR:  Times do not match!  Aborting...'
        WRITE(6,*)'mtimesec = ',mtimesec,' timesec = ',timesec
        STOP
      END IF
    END DO

!-----------------------------------------------------------------------
!
!   Copy the bias and rms scores to the "master" arrays, and then move
!   to the next file.
! 
!-----------------------------------------------------------------------

    mbias_p(:,:,:,i) = bias_p
    mrms_p(:,:,:,i) = rms_p
    mbias_sfc(:,:,i) = bias_sfc
    mrms_sfc(:,:,i) = rms_sfc
    mcounter(i) = counter
    mcounter_p(:,:,:,i) = counter_p
    mcounter_sfc(:,:,i) = counter_sfc
    mdiff_p(:,:,:,:,:,i) = diff_p
    mdiff_sfc(:,:,:,:,i) = diff_sfc

  END DO

!-----------------------------------------------------------------------
!
! Now calculate the average bias and rmse scores.
! 
!-----------------------------------------------------------------------

  DO i = 1, numinfiles
    IF (mcounter(i).NE.missing .AND. mcounter(i).NE.0) THEN
      avgcounter = avgcounter + mcounter(i)
      avgcounterfiles = avgcounterfiles + 1
    END IF
  END DO

  DO i = 1,nlevel
    DO j = 1,maxntime
      DO k = 1,nvar_p_stats
        DO l = 1,numinfiles
          IF (mbias_p(i,j,k,l).NE.missing .AND. &
              mrms_p(i,j,k,l).NE.missing) THEN
            avgbias_p(i,j,k) = avgbias_p(i,j,k) + mbias_p(i,j,k,l)
            avgbiasfiles_p(i,j,k) = avgbiasfiles_p(i,j,k) + 1
            avgrms_p(i,j,k) = avgrms_p(i,j,k) + mrms_p(i,j,k,l)
            avgrmsfiles_p(i,j,k) = avgrmsfiles_p(i,j,k) + 1
            avgcounter_p(i,j,k) = avgcounter_p(i,j,k) + &
              mcounter_p(i,j,k,l)
            avgcounterfiles_p(i,j,k) = avgcounterfiles_p(i,j,k) + 1
          END IF          
        END DO
      END DO
    END DO
  END DO

  DO i = 1,maxntime
    DO j = 1,nvar_sfc_stats
      DO k = 1,numinfiles
        IF (mbias_sfc(i,j,k).NE.missing .AND. &
            mrms_sfc(i,j,k).NE.missing) THEN
          avgbias_sfc(i,j) = avgbias_sfc(i,j) + mbias_sfc(i,j,k)
          avgbiasfiles_sfc(i,j) = avgbiasfiles_sfc(i,j) + 1
          avgrms_sfc(i,j) = avgrms_sfc(i,j) + mrms_sfc(i,j,k)
          avgrmsfiles_sfc(i,j) = avgrmsfiles_sfc(i,j) + 1
          avgcounter_sfc(i,j) = avgcounter_sfc(i,j) + &
            mcounter_sfc(i,j,k)
          avgcounterfiles_sfc(i,j) = avgcounterfiles_sfc(i,j) + 1
        END IF
      END DO
    END DO
  END DO

  DO i = 1,nlevel
    DO j = 1,maxntime
      DO k = 1,nvar_p_stats
        IF (avgbiasfiles_p(i,j,k).GT.0) THEN
          avgbias_p(i,j,k) = avgbias_p(i,j,k)/REAL(avgbiasfiles_p(i,j,k))
        ELSE
          avgbias_p(i,j,k) = missing
        END IF

        IF (avgrmsfiles_p(i,j,k).GT.0) THEN
          avgrms_p(i,j,k) = avgrms_p(i,j,k)/REAL(avgrmsfiles_p(i,j,k))
        ELSE
          avgrms_p(i,j,k) = missing
        END IF

        IF (avgcounterfiles_p(i,j,k).GT.0) THEN
          avgcounter_p(i,j,k) = &
            avgcounter_p(i,j,k)/REAL(avgcounterfiles_p(i,j,k))
        ELSE
          avgcounter_p(i,j,k) = 0
        END IF

      END DO
    END DO
  END DO

  DO i = 1,maxntime
    DO j = 1,nvar_sfc_stats

      IF (avgbiasfiles_sfc(i,j).GT.0) THEN
        avgbias_sfc(i,j) = avgbias_sfc(i,j)/REAL(avgbiasfiles_sfc(i,j))
      ELSE
        avgbias_sfc(i,j) = missing
      END IF

      IF (avgrmsfiles_sfc(i,j).GT.0) THEN
        avgrms_sfc(i,j) = avgrms_sfc(i,j)/REAL(avgrmsfiles_sfc(i,j))
      ELSE
        avgrms_sfc(i,j) = missing
      END IF

      IF (avgcounterfiles_sfc(i,j).GT.0) THEN
        avgcounter_sfc(i,j) = &
          avgcounter_sfc(i,j)/REAL(avgcounterfiles_sfc(i,j))
      ELSE
        avgcounter_sfc(i,j) = 0
      END IF

    END DO
  END DO

  IF (avgcounterfiles.GT.0) THEN
    avgcounter = avgcounter/avgcounterfiles
  ELSE
    avgcounter = 0
  END IF

!-----------------------------------------------------------------------
!
! Now write the average statstics to an output file.
! 
!-----------------------------------------------------------------------

  PRINT *, 'Writing text output to ', TRIM(txtout)
  OPEN(UNIT=13,FILE=txtout,STATUS="REPLACE")
      
  590 FORMAT (1x,a,i2.2,a)
  WRITE(13,590) 'Average statistics for ',numinfiles,' forecast runs.'
  WRITE(13,610) 'Forecast model/grid:', TRIM(fmodel), TRIM(fgrid)
  WRITE(13,626) fdx/1000.0

  WRITE(13,610) 'Analysis model/grid:', TRIM(amodel), TRIM(agrid)
  WRITE(13,626)  adx/1000.0
  WRITE(13,*)
  WRITE(13,*) 'Verification Region:' 
  WRITE(13,625) NINT((vnx+1)*vdx/1000.0), NINT((vny+1)*vdx/1000.0), &
                vdx/1000.0, vnx,vny
  WRITE(13,650) 'NW/NE', &
      corner_lat(1,2),corner_lon(1,2),corner_lat(2,2),corner_lon(2,2)
  WRITE(13,650) 'SW/SE', &
      corner_lat(1,1),corner_lon(1,1),corner_lat(2,1),corner_lon(2,1)
  WRITE(13,*) 'Average number of grid points in verif region: ', &
              avgcounter
  WRITE(13,*)

  DO i = 1,numinfiles
    WRITE(13,600) 'Forecast initialized ',mfdate(1,i),'-',mfdate(2,i),&  
      '-',mfdate(3,i),'.',mfdate(4,i),':',mfdate(5,i),':',mfdate(6,i)
  END DO
  
  DO l = 1,nvar_sfc_stats
    WRITE(13,'(A)') TRIM(line72)
    WRITE(13,810) TRIM(varname_sfc(l)), TRIM(varunit_sfc(l))
    WRITE(13,*) 'Time(hr)  Grid Points      Bias  RMS Error',&
                '  Fcsts'
    DO k = 1,maxntime
      WRITE(13,700) timesec(k)/3600, avgcounter_sfc(k,l), &
          avgbias_sfc(k,l), avgrms_sfc(k,l), avgbiasfiles_sfc(k,l)
    END DO
  END DO
  
  DO l = 1,nvar_p_stats
  DO m = 1,nlevel
    WRITE(13,'(A)') TRIM(line72)
    WRITE(13,810) TRIM(varname_p(l)), TRIM(varunit_p(l))
    WRITE(13,*) 'Level = ', pressure(m) / 100.0, ' mb'
    WRITE(13,*) 'Time(hr)  Grid Points      Bias  RMS Error',&
                '  Fcsts'
    DO k = 1,maxntime
      WRITE(13,700) timesec(k)/3600, avgcounter_p(m,k,l), &
          avgbias_p(m,k,l), avgrms_p(m,k,l),avgbiasfiles_p(m,k,l)
    END DO
  END DO
  END DO

  CALL wrt_verif_stats_diff (hdfout, 0, hdiffile, mfdate(:,1), &
                     mamodel,magrid,madx,anx,any, &
                     mfmodel,mfgrid,mfdx, &
                     mvdx, mvdy, mvnx, mvny, mvmapproj, &
                     mvtrulat1, mvtrulat2, mvtrulon, mvsclfct, &
                     mvctrlat, mvctrlon, &
                     mcorner_lat, mcorner_lon, &
                     nlevel, maxntime, nvar_p_stats, nvar_sfc_stats, &  
                     timesec, pressure, missing, avgcounter, &
                     avgcounter_p, avgbias_p, avgrms_p, diff_p, &
                     avgcounter_sfc, avgbias_sfc, avgrms_sfc, diff_sfc, &
                     varid_p, varname_p, varunit_p, &
                     varid_sfc, varname_sfc, varunit_sfc)

  600 FORMAT (1x,a21,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)
  610 FORMAT (1X,A,1X,A,' / ',A)
  625 FORMAT (' Domain:',I5, ' x', I5, ' km, dx=', F5.1, &
              ' km (',I4,' x',I4,')' )
  626 FORMAT (' dx=', F5.1, ' km')
  650 FORMAT (1X,A,' corner lat/lon:', 2F10.3, 2X, 2F10.3)
  700 FORMAT (1x,i8,2x,i11,2x,f8.2,3x,f8.2,3x,I4)
  800 FORMAT (1x,i8,2x,i11,4x,f8.2,4x,f8.2,10x,f8.2)
  810 FORMAT ( 1X,A," (",A,")" )

END PROGRAM avgverifgrid
