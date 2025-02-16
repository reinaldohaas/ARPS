!########################################################################
!########################################################################
!######                                                            ######
!######                    PROGRAM QPFSTATS                        ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

PROGRAM QPFSTATS

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Reads in two HDF files created by INTQPF (a forecast and a
! verification) and a binary bit map created by QPFMASK.  The precip.
! fields in the HDF files extracted, bias and equitable threat score
! are calculated, and the statistics are output.
!
! AUTHOR:  Eric Kemp, March 2000
!
! MODIFICATION HISTORY:
! Eric Kemp, 31 March 2000.
! Added lat/lon coordinates of four corners for NCL plotting.
!
!-----------------------------------------------------------------------
! 
! Variable declarations
! 
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! 
! Forecast variables 
! 
!-----------------------------------------------------------------------

  INTEGER :: fdate(6)
  INTEGER,ALLOCATABLE :: ftimesec(:)
  
  INTEGER :: fmapproj
  REAL :: fscale
  REAL :: ftrulat1, ftrulat2, ftrulon
  REAL :: fdx,fdy
                       
  CHARACTER*8 :: fmodel
  CHARACTER*4 :: fgrid
  CHARACTER*4,ALLOCATABLE,DIMENSION(:) :: fid_p,fid_sfc
  CHARACTER*32,ALLOCATABLE,DIMENSION(:) :: fname_p,funit_p,fname_sfc, &
                                           funit_sfc

  REAL,ALLOCATABLE, DIMENSION(:,:,:,:,:) :: fvar_p
  REAL,ALLOCATABLE, DIMENSION(:,:,:,:) :: fvar_sfc
  REAL,ALLOCATABLE, DIMENSION(:) :: fpressure
  REAL,ALLOCATABLE, DIMENSION(:,:) :: flat2d,flon2d
  REAL :: fscswlat,fscswlon
  REAL :: fmapswlat,fmapswlon,fmapnwlat,fmapnwlon, &
          fmapselat,fmapselon,fmapnelat,fmapnelon

  INTEGER :: fflag_p,fflag_sfc

  REAL,ALLOCATABLE, DIMENSION(:,:,:) :: fprecip

  INTEGER :: fnx,fny,fntime
  INTEGER :: fnpreslevel,fnvar_p,fnvar_sfc

!-----------------------------------------------------------------------
! 
! Verification variables 
! 
!-----------------------------------------------------------------------

  INTEGER :: vdate(6)
  INTEGER,ALLOCATABLE :: vtimesec(:)
  
  INTEGER :: vmapproj
  REAL :: vscale
  REAL :: vtrulat1, vtrulat2, vtrulon
  REAL :: vdx,vdy
                       
  CHARACTER*8 :: vmodel
  CHARACTER*4 :: vgrid
  CHARACTER*4,ALLOCATABLE,DIMENSION(:) :: vid_p,vid_sfc
  CHARACTER*32,ALLOCATABLE,DIMENSION(:) :: vname_p,vunit_p,vname_sfc, &
                                           vunit_sfc

  REAL,ALLOCATABLE, DIMENSION(:,:,:,:,:) :: vvar_p
  REAL,ALLOCATABLE, DIMENSION(:,:,:,:) :: vvar_sfc
  REAL,ALLOCATABLE, DIMENSION(:) :: vpressure
  REAL,ALLOCATABLE, DIMENSION(:,:) :: vlat2d,vlon2d
  REAL :: vscswlat,vscswlon
  REAL :: vmapswlat,vmapswlon,vmapnwlat,vmapnwlon, &
          vmapselat,vmapselon,vmapnelat,vmapnelon
  
  INTEGER :: vflag_p,vflag_sfc

  REAL,ALLOCATABLE, DIMENSION(:,:,:) :: vprecip
  INTEGER :: vnx,vny,vntime
  INTEGER :: vnpreslevel,vnvar_p,vnvar_sfc
  
!-----------------------------------------------------------------------
! 
! Bitmap variables 
! 
!-----------------------------------------------------------------------

  INTEGER :: mdate(6)
  INTEGER,ALLOCATABLE :: mtimesec(:)
  
  INTEGER :: mmapproj
  REAL :: mscale
  REAL :: mtrulat1, mtrulat2, mtrulon
  REAL :: mdx,mdy
                       
  CHARACTER*8 :: mmodel
  CHARACTER*4 :: mgrid
  CHARACTER*4,ALLOCATABLE,DIMENSION(:) :: mid_p,mid_sfc
  CHARACTER*32,ALLOCATABLE,DIMENSION(:) :: mname_p,munit_p,mname_sfc, &
                                           munit_sfc

  REAL,ALLOCATABLE, DIMENSION(:,:,:,:,:) :: mvar_p
  REAL,ALLOCATABLE, DIMENSION(:,:,:,:) :: mvar_sfc
  REAL,ALLOCATABLE, DIMENSION(:) :: mpressure
  REAL,ALLOCATABLE, DIMENSION(:,:) :: mlat2d,mlon2d
  REAL :: mscswlat,mscswlon
  REAL :: mmapswlat,mmapswlon,mmapnwlat,mmapnwlon, &
          mmapselat,mmapselon,mmapnelat,mmapnelon

  INTEGER :: mflag_p,mflag_sfc

  INTEGER,ALLOCATABLE, DIMENSION(:,:) :: mask
  INTEGER :: mnx,mny,mntime
  INTEGER :: mnpreslevel,mnvar_p,mnvar_sfc

  INTEGER :: found_mask

!-----------------------------------------------------------------------
! 
! Other variables
!
! NOTE:  CTA = array of "yes" forecast, "yes" observation grid points
!        CTB = array of "yes" forecast, "no" observation grid points
!        CTC = array of "no" forecast, "yes" observation grid points
!        CTD = array of "no" forecast, "no" observation grid points
! 
! Reference:  Wilks, D. S., 1995:  Statistical Methods in the 
!             Atmospheric Sciences.  Academic Press, 476pp.
!
!-----------------------------------------------------------------------
  
  REAL,ALLOCATABLE :: ets(:,:),bias(:,:)
  INTEGER,ALLOCATABLE :: CTA(:,:),CTB(:,:),CTC(:,:),CTD(:,:)

  INTEGER :: status,i,j,k,l
  REAL :: factor
  REAL,PARAMETER :: missing = -9999.

  INTEGER,PARAMETER :: nstat = 2
  CHARACTER*4,ALLOCATABLE,DIMENSION(:) :: ctstatsid
  CHARACTER*32,ALLOCATABLE,DIMENSION(:) :: ctstatsname
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ctstats
  CHARACTER*32 :: threshunit

  CHARACTER*72,PARAMETER :: line72 = &
 '------------------------------------------------------------------------'

  REAL :: vtrulat(2)
  REAL,ALLOCATABLE :: vxsc(:),vysc(:)
  REAL :: vx0,vy0

!-----------------------------------------------------------------------
!
! Namelists
!
!-----------------------------------------------------------------------
  
  CHARACTER*256 :: finfilename
  NAMELIST /hdf_finput/ finfilename

  CHARACTER*256 :: vinfilename
  NAMELIST /hdf_vinput/ vinfilename

  CHARACTER*256 :: minfilename
  NAMELIST /bitmap/ minfilename

  INTEGER :: precipunit
  INTEGER :: numthresh
  INTEGER,PARAMETER :: maxnumthresh = 10
  REAL :: thresh(maxnumthresh)
  NAMELIST /threshold/ precipunit,numthresh,thresh

  CHARACTER*256 :: hdfoutfilename,ascoutfilename
  NAMELIST /output/ hdfoutfilename,ascoutfilename

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  fmodel=''
  vmodel=''

  WRITE(6,'(//5x,a)')                                                   &
 &'###################################################################'
  WRITE(6,'(5x,a,/5x,a)')                                               &
 &'#                                                                 #',&
 &'# Welcome to QPFSTATS, a program that reads in two precip. HDF    #'
  WRITE(6,'(5x,a)')                                                     &
 &'# files from INTQPF and a binary bitmap file from QPFMASK, and    #'
  WRITE(6,'(5x,a)')                                                     &
 &'# calculates equitable threat score and bias statistics.          #', &
 &'#                                                                 #',&
 &'###################################################################'
  WRITE(6,*)

!-----------------------------------------------------------------------
!
! Read in namelist
!
!----------------------------------------------------------------------- 
  
  PRINT *, 'Reading NAMELIST hdf_finput'
  READ(5,hdf_finput)
  PRINT *, 'Reading NAMELIST hdf_vinput'
  READ(5,hdf_vinput)
  PRINT *, 'Reading NAMELIST bitmap'
  READ(5,bitmap)
  PRINT *, 'Reading NAMELIST threshold'
  READ(5,threshold)
  PRINT *, 'Reading NAMELIST output'
  READ(5,output)

  WRITE(6,*)'EMK: thresh = ',thresh

  IF (numthresh > maxnumthresh) THEN
    WRITE(6,*)'ERROR:  numthresh > maxnumthresh'
    WRITE(6,*)'Change maxnumthresh and recompile.'
    WRITE(6,*)'numthresh = ',numthresh,' maxnumthresh = ',   &
              maxnumthresh
    WRITE(6,*)'Aborting...'
    STOP 
  END IF

  IF (precipunit /= 1 .AND. precipunit /= 2) THEN
    WRITE(6,*)'ERROR:  Invalid precipunit value'
    WRITE(6,*)'precipunit = ',precipunit
    WRITE(6,*)'Aborting...'    
    STOP
  ENDIF

!----------------------------------------------------------------------- 
!
! Read forecast precipitation data.
!
!-----------------------------------------------------------------------

  CALL rd_verif_dims(fnx,fny,fnpreslevel,fntime,fnvar_p,fnvar_sfc, &
                     finfilename)

  ALLOCATE (ftimesec(fntime), fpressure(fnpreslevel),        &
            fvar_p(fnx,fny,fnpreslevel,fntime,fnvar_p), &
            fvar_sfc(fnx,fny,fntime,fnvar_sfc), &
            flat2d(fnx,fny),flon2d(fnx,fny),  &
            fprecip(fnx,fny,fntime), &
            fid_p(fnvar_p),fname_p(fnvar_p),funit_p(fnvar_p), &
            fid_sfc(fnvar_sfc),fname_sfc(fnvar_sfc),          &
            funit_sfc(fnvar_sfc),  &
            STAT=status)

  IF (status /= 0) CALL alloc_fail (status, 'f1')

  CALL rdverif ( fnx,fny,fnpreslevel,fntime,fnvar_p,fnvar_sfc,      &
                 finfilename,fmodel,fgrid,fdate,ftimesec,fpressure, &
                 fmapproj,fscale,ftrulat1,ftrulat2,ftrulon,fdx,fdy, &
                 fscswlat,fscswlon, &
                 fmapswlat,fmapswlon,fmapnwlat,fmapnwlon, &
                 fmapselat,fmapselon,fmapnelat,fmapnelon, &
                 fflag_p,fvar_p, fid_p, fname_p, funit_p, &
                 fflag_sfc,fvar_sfc, fid_sfc, fname_sfc, funit_sfc )

  IF (fflag_sfc == 0) THEN 
    WRITE(6,*)'ERROR:  Could not find surface data in HDF file.'
    WRITE(6,*)'Aborting...'
    STOP
  ENDIF

!----------------------------------------------------------------------- 
!
! Read verification precipitation data.
!
!-----------------------------------------------------------------------

  CALL rd_verif_dims(vnx,vny,vnpreslevel,vntime,vnvar_p,vnvar_sfc, &
                     vinfilename)

  ALLOCATE (vtimesec(vntime), vpressure(vnpreslevel),        &
            vvar_p(vnx,vny,vnpreslevel,vntime,vnvar_p), &
            vvar_sfc(vnx,vny,vntime,vnvar_sfc), &
            vlat2d(vnx,vny),vlon2d(vnx,vny),  &
            vprecip(vnx,vny,vntime),  &
            vid_p(vnvar_p),vname_p(vnvar_p),vunit_p(vnvar_p), &
            vid_sfc(vnvar_sfc),vname_sfc(vnvar_sfc),          &
            vunit_sfc(vnvar_sfc),  &
            STAT=status)

  IF (status /= 0) CALL alloc_fail (status, 'f1')

  CALL rdverif ( vnx,vny,vnpreslevel,vntime,vnvar_p,vnvar_sfc,      &
                 vinfilename,vmodel,vgrid,vdate,vtimesec,vpressure, &
                 vmapproj,vscale,vtrulat1,vtrulat2,vtrulon,vdx,vdy, &
                 vscswlat,vscswlon, &
                 vmapswlat,vmapswlon,vmapnwlat,vmapnwlon, &
                 vmapselat,vmapselon,vmapnelat,vmapnelon, &
                 vflag_p,vvar_p, vid_p, vname_p, vunit_p, &
                 vflag_sfc,vvar_sfc, vid_sfc, vname_sfc, vunit_sfc )

  IF (vflag_sfc == 0) THEN 
    WRITE(6,*)'ERROR:  Could not find surface data in HDF file.'
    WRITE(6,*)'Aborting...'
    STOP
  ENDIF

!----------------------------------------------------------------------- 
!
! Check dimensions, map projection, etc. of forecast and verification.
!
!-----------------------------------------------------------------------

  IF (fnx /= vnx  .OR. fny /= vny .OR. fntime /= vntime .OR.        &
      fscale /= vscale .OR.                     &
      ftrulat1 /= vtrulat1 .OR. ftrulat2 /= vtrulat2 .OR.           &
      ftrulon /= vtrulon .OR. fdx /= vdx .OR. fdy /= vdy .OR. &
      fscswlat /= vscswlat .OR. fscswlon /= vscswlon) THEN

    WRITE(6,*)'ERROR:  Grid information for forecast and verification'
    WRITE(6,*)'data do not match!'
    WRITE(6,*)'fnx = ',fnx,' vnx = ',vnx    
    WRITE(6,*)'fny = ',fny,' vny = ',vny    
    WRITE(6,*)'fscale = ',fscale,' vscale = ',vscale
    WRITE(6,*)'ftrulat1 = ',ftrulat1,' vtrulat1 = ',vtrulat1
    WRITE(6,*)'ftrulat2 = ',ftrulat2,' vtrulat2 = ',vtrulat2
    WRITE(6,*)'ftrulon = ',ftrulon,' vtrulon = ',vtrulon
    WRITE(6,*)'fdx = ',fdx,' vdx = ',vdx
    WRITE(6,*)'fdy = ',fdy,' vdy = ',vdy
    WRITE(6,*)'fscswlat = ',fscswlat,' vscswlat = ',vscswlat
    WRITE(6,*)'fscswlon = ',fscswlon,' vscswlon = ',vscswlon
    WRITE(6,*)
    WRITE(6,*)'Aborting...'
    STOP 
  END IF

  DO i=1,6
    IF (fdate(i) /= vdate(i)) THEN
      PRINT *, 'FATAL: Forecast and analysis start dates do not match'
      PRINT *, 'fdate= ', fdate
      PRINT *, 'vdate= ', vdate
    STOP
    END IF
  END DO

!----------------------------------------------------------------------- 
!
! Read bitmap file.
!
!-----------------------------------------------------------------------

  CALL rd_verif_dims(mnx,mny,mnpreslevel,mntime,mnvar_p,mnvar_sfc, &
                     minfilename)

  ALLOCATE (mtimesec(mntime), mpressure(mnpreslevel),        &
            mvar_p(mnx,mny,mnpreslevel,mntime,mnvar_p), &
            mvar_sfc(mnx,mny,mntime,mnvar_sfc), &
            mlat2d(mnx,mny),mlon2d(mnx,mny),  &
            mid_p(mnvar_p),mname_p(mnvar_p),munit_p(mnvar_p), &
            mid_sfc(mnvar_sfc),mname_sfc(mnvar_sfc),          &
            munit_sfc(mnvar_sfc),  &
            STAT=status)

  IF (status /= 0) CALL alloc_fail (status, 'f1')

  CALL rdverif ( mnx,mny,mnpreslevel,mntime,mnvar_p,mnvar_sfc,      &
                 minfilename,mmodel,mgrid,mdate,mtimesec,mpressure, &
                 mmapproj,mscale,mtrulat1,mtrulat2,mtrulon,mdx,mdy, &
                 mscswlat,mscswlon, &
                 mmapswlat,mmapswlon,mmapnwlat,mmapnwlon, &
                 mmapselat,mmapselon,mmapnelat,mmapnelon, &
                 mflag_p,mvar_p, mid_p, mname_p, munit_p, &
                 mflag_sfc,mvar_sfc, mid_sfc, mname_sfc, munit_sfc )

  IF (mflag_sfc == 0) THEN 
    WRITE(6,*)'ERROR:  Could not find bitmap data in HDF file.'
    WRITE(6,*)'Aborting...'
    STOP
  ENDIF

!----------------------------------------------------------------------- 
!
! Check dimensions, map projection, etc. of bitmap and verification
!
!-----------------------------------------------------------------------

  IF (mnx /= vnx  .OR. mny /= vny .OR. mntime /= vntime .OR.        &
      mscale /= vscale .OR.                     &
      mtrulat1 /= vtrulat1 .OR. mtrulat2 /= vtrulat2 .OR.           &
      mtrulon /= vtrulon .OR. mdx /= vdx .OR. mdy /= vdy .OR. &
      mscswlat /= vscswlat .OR. mscswlon /= vscswlon) THEN 

    WRITE(6,*)'ERROR:  Grid information for bitmap and verification'
    WRITE(6,*)'data do not match!'
    WRITE(6,*)'mnx = ',mnx,' vnx = ',vnx    
    WRITE(6,*)'mny = ',mny,' vny = ',vny    
    WRITE(6,*)'mscale = ',mscale,' vscale = ',vscale
    WRITE(6,*)'mtrulat1 = ',mtrulat1,' vtrulat1 = ',vtrulat1
    WRITE(6,*)'mtrulat2 = ',mtrulat2,' vtrulat2 = ',vtrulat2
    WRITE(6,*)'mtrulon = ',mtrulon,' vtrulon = ',vtrulon
    WRITE(6,*)'mdx = ',mdx,' vdx = ',vdx
    WRITE(6,*)'mdy = ',mdy,' vdy = ',vdy
    WRITE(6,*)'mscswlat = ',mscswlat,' vscswlat = ',vscswlat
    WRITE(6,*)'mscswlon = ',mscswlon,' vscswlon = ',vscswlon
    WRITE(6,*)
    WRITE(6,*)'Aborting...'
    STOP 
  END IF

!----------------------------------------------------------------------- 
!
! Dump HDF data into precip arrays, changing units as necessary.
!
!-----------------------------------------------------------------------
  
  DO i = 1,fnvar_sfc
    IF (fid_sfc(i) == 'APCP') THEN
      IF (funit_sfc(i) == 'mm' .AND. precipunit == 1) THEN
        fprecip = fvar_sfc(:,:,:,i)      
      ELSE IF (funit_sfc(i) == 'mm' .AND. precipunit == 2) THEN
        fprecip = fvar_sfc(:,:,:,i)
        factor = 39.37/1000.
        CALL cvtpcpunit(fnx,fny,fntime,fprecip,factor,missing)
      ELSE IF (funit_sfc(i) == 'in' .AND. precipunit == 2) THEN
        fprecip = fvar_sfc(:,:,:,i)
      ELSE IF (funit_sfc(i) == 'in' .AND. precipunit == 1) THEN
        fprecip = fvar_sfc(:,:,:,i)
        factor = 1000./39.37
        CALL cvtpcpunit(fnx,fny,fntime,fprecip,factor,missing)
      ELSE 
        WRITE(6,*)'ERROR:  funit_sfc = ',funit_sfc
        STOP
      ENDIF      
    ENDIF
  END DO

  DO i = 1,vnvar_sfc
    IF (vid_sfc(i) == 'APCP') THEN
      IF (vunit_sfc(i) == 'mm' .AND. precipunit == 1) THEN
        vprecip = vvar_sfc(:,:,:,i)
      ELSE IF (vunit_sfc(i) == 'mm' .AND. precipunit == 2) THEN
        vprecip = vvar_sfc(:,:,:,i)
        factor = 39.37/1000.
        CALL cvtpcpunit(vnx,vny,vntime,vprecip,factor,missing)
      ELSE IF (vunit_sfc(i) == 'in' .AND. precipunit == 2) THEN
        vprecip = vvar_sfc(:,:,:,i)
      ELSE IF (vunit_sfc(i) == 'in' .AND. precipunit == 1) THEN
        vprecip = vvar_sfc(:,:,:,i)
        factor = 1000./39.37
        CALL cvtpcpunit(vnx,vny,vntime,vprecip,factor,missing)
      ELSE 
        WRITE(6,*)'ERROR:  vunit_sfc = ',vunit_sfc
        STOP
      ENDIF      
    ENDIF
  END DO

!----------------------------------------------------------------------- 
!
! Dump HDF bitmap data into mask array.
!
!-----------------------------------------------------------------------

  ALLOCATE (mask(mnx,mny),STAT=status)
  IF (status /= 0) CALL alloc_fail (status, 'f1')

  mask = missing

  found_mask = 0

  DO i = 1,mnvar_sfc
    IF (mid_sfc(i) == 'MASK') THEN
      mask = INT(mvar_sfc(:,:,1,i))
      found_mask = 1
    END IF
  END DO

  IF (found_mask == 0) THEN
    WRITE(6,*)'ERROR:  Could not find bitmap.  Aborting...'
    STOP
  END IF

!----------------------------------------------------------------------- 
!
! Calculate equitable threat score and bias
!
!-----------------------------------------------------------------------

  ALLOCATE (bias(vntime,numthresh),ets(vntime,numthresh),  &
            CTA(vntime,numthresh),CTB(vntime,numthresh), &
            CTC(vntime,numthresh),CTD(vntime,numthresh), &
            STAT=status)
  IF (status /= 0) CALL alloc_fail (status, 'f1')

  CALL calcetsb(vnx,vny,vntime,vprecip,fprecip,numthresh,thresh, &
                mask,bias,ets,CTA,CTB,CTC,CTD)

  100 FORMAT(1x,a,f4.2,a)
  110 FORMAT(1x,A,i5,A,i5,A,i5)

!----------------------------------------------------------------------- 
!
! Dump statistical data to a new HDF file.
!
!-----------------------------------------------------------------------

  ALLOCATE(ctstats(vntime,numthresh,nstat),ctstatsid(nstat),  &
           ctstatsname(nstat),STAT=status)
  IF (status /= 0) CALL alloc_fail (status, 'f1')
 
  ctstatsid(1) = 'BIAS'
  ctstatsid(2) = 'ETS '
  ctstatsname(1) = 'Bias'
  ctstatsname(2) = 'Equitable Threat Score'

  IF (precipunit == 1) THEN
    threshunit = 'mm'
  ELSE IF (precipunit == 2) THEN
    threshunit = 'in'
  ELSE
    WRITE(6,*)'ERROR:  Invalid precipunit value'
    WRITE(6,*)'precipunit = ',precipunit
    WRITE(6,*)'Aborting...'
    STOP
  END IF

  DO i = 1,nstat
    IF (i == 1) THEN
      ctstats(:,:,i) = bias
    ELSE IF (i == 2) THEN
      ctstats(:,:,i) = ets
    ENDIF
  END DO
  CALL wrtqpfstats(vnx,vny,vntime,numthresh,nstat,missing, &
                   hdfoutfilename,vmodel,vgrid, &
                   vdate,vtimesec,vmapproj,vscale,vtrulat1,vtrulat2, &
                   vtrulon,vdx,vdy,vscswlat,vscswlon,&
                   vmapswlat,vmapswlon,vmapnwlat,vmapnwlon,             &
                   vmapselat,vmapselon,vmapnelat,vmapnelon,             &
                   thresh,threshunit,  &
                   ctstatsid,ctstatsname,ctstats,                    &
                   CTA,CTB,CTC,CTD)

!----------------------------------------------------------------------- 
!
! Output ASCII data.
!
!-----------------------------------------------------------------------

  ALLOCATE(vxsc(vnx),vysc(vny),STAT=status)
  IF (status /= 0) CALL alloc_fail (status, 'f1')

  vtrulat(1) = vtrulat1 
  vtrulat(2) = vtrulat2
 
  CALL setmapr(vmapproj,vscale,vtrulat,vtrulon)
  CALL lltoxy(1,1,vscswlat,vscswlon,vx0,vy0)
    
  DO i=1,vnx
    vxsc(i)=vx0+(i-1)*vdx
  END DO
  
  DO j=1,vny
    vysc(j)=vy0+(j-1)*vdy
  END DO
  
  CALL xytoll(vnx,vny,vxsc,vysc,vlat2d,vlon2d)


  WRITE(6,*)'Writing text output to ',TRIM(ascoutfilename)
  OPEN(UNIT=13,FILE=ascoutfilename,STATUS="REPLACE")
  
  WRITE(13,600) 'Forecast initialized ',fdate(1),'-',fdate(2),        &
            '-',fdate(3),'.',fdate(4),':',fdate(5),':',fdate(6)
  WRITE(13,610) 'Forecast precipitation data:', TRIM(fmodel)
  WRITE(13,610) 'Analysis precipitation data:', TRIM(vmodel)

  WRITE(13,*)
  WRITE(13,*) 'Verification Region:' 
  WRITE(13,625) NINT((vnx)*vdx/1000.0), NINT((vny)*vdx/1000.0), &
                vdx/1000.0, vnx,vny
  WRITE(13,650) 'NW/NE', &
      vlat2d(1,vny),vlon2d(1,vny),vlat2d(vnx,vny),vlon2d(vnx,vny)
  WRITE(13,650) 'SW/SE', &
      vlat2d(1,1),vlon2d(1,1), vlat2d(vnx,1),vlon2d(vnx,1)
  WRITE(13,*)'Used bit map from file ',TRIM(minfilename)
  WRITE(13,*)
  
  DO l = 1, numthresh
    WRITE(13,*) line72
    WRITE(13,810) 'Precipitation Threshold:  ',thresh(l), TRIM(threshunit)
    WRITE(13,*) 'Time    Grid                 False           '
    WRITE(13,*) '(hr)  Points  Hits  Misses  Alarms      Bias', &
    '        ETS'
    DO k = 1,vntime
      WRITE(13,700) vtimesec(k)/3600, &
                   CTA(k,l)+CTB(k,l)+CTC(k,l)+CTD(k,l), &
                   CTA(k,l),CTC(k,l),CTB(k,l), &
                   bias(k,l),ETS(k,l)
    END DO
  END DO

  600 FORMAT (/1x,a21,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)
  610 FORMAT (1X,A,1X,A)
  625 FORMAT (' Domain:',I5, ' x', I5, ' km, dx=', F5.1, &
              ' km (',I4,' x',I4,' )' )
  650 FORMAT (1X,A,' corner lat/lon:', 2F10.3, 2X, 2F10.3)
  700 FORMAT (3x,i2,4x,i4,2x,i4,4x,i4,4x,i4,2x,f8.2,3x,f8.2)
  810 FORMAT ( 1X,a,f4.2," (",A,")" )
  
!----------------------------------------------------------------------- 
!
! The end.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'Program QPFSTATS successfully completed.'
 
END PROGRAM QPFSTATS


!########################################################################
!########################################################################
!######                                                            ######
!######                  SUBROUTINE CALCETSB                       ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE calcetsb(nx,ny,ntime,vprecip,fprecip,numthresh,thresh, &
                    mask,bias,ets,cta,ctb,ctc,ctd)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculates equitable threat score and bias scores for precipitation.
!
! AUTHOR:  Eric Kemp, March 2000
!
!-----------------------------------------------------------------------
! 
! Variable declarations
! 
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,ntime,numthresh
  REAL, INTENT(IN) :: vprecip(nx,ny,ntime), fprecip(nx,ny,ntime)
  REAL, INTENT(OUT) :: ETS(ntime,numthresh),Bias(ntime,numthresh)
  REAL, INTENT(IN) :: thresh(numthresh)
  INTEGER, INTENT(IN) :: mask(nx,ny)
  INTEGER, INTENT(INOUT), DIMENSION(ntime,numthresh) :: CTA,CTB,CTC,CTD

!-----------------------------------------------------------------------
! 
! Variables for statistical calculations 
!
! LA = Number of "yes" forecast, "yes" observation grid points.
! LB = Number of "yes" forecast, "no" observation grid points.
! LC = Number of "no" forecast, "yes" observation grid points.
! LD = Number of "no" forecast, "no" observation grid points.
! R = Number of random correct forecasts expected due to chance within
!     the total number of model and observation pairs.
! ETS = Equitable Threat Score.
!
! NOTE:  R = (LA + LB)*(LA + LC)/(LA + LB + LC + LD)
!        ETS = (LA - R)/(LA + LB + LC - R)
!        Bias = (LA + LB)/(LA + LC)
!
!-----------------------------------------------------------------------

  REAL :: R
  INTEGER :: LA, LB, LC, LD

!-----------------------------------------------------------------------
! 
! Other variables
! 
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,l

  REAL, PARAMETER :: missing = -9999.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Bias = missing
  ETS = missing
  CTA = 0
  CTB = 0
  CTC = 0
  CTD = 0

  DO l = 1,numthresh
    DO k = 1,ntime
      LA = 0
      LB = 0
      LC = 0
      LD = 0
      DO j = 1,ny
        DO i = 1,nx
          IF (mask(i,j) == 1 .AND. fprecip(i,j,k) /= missing .AND.  &
              vprecip(i,j,k) /= missing ) THEN
            IF (fprecip(i,j,k) >= thresh(l) .AND.  &
                vprecip(i,j,k) >= thresh(l)) THEN
              LA = LA + 1
            ELSE IF (fprecip(i,j,k) >= thresh(l) .AND.  &
                vprecip(i,j,k) < thresh(l)) THEN
              LB = LB + 1
            ELSE IF (fprecip(i,j,k) < thresh(l) .AND.  &
                vprecip(i,j,k) >= thresh(l)) THEN
              LC = LC + 1
            ELSE 
              LD = LD + 1
            ENDIF
          ENDIF
        END DO          
      END DO          

      IF ((LA + LB + LC + LD) == 0) THEN
        R = missing
      ELSE
        R = REAL(LA + LB)*REAL(LA + LC)/REAL(LA + LB + LC + LD)
      ENDIF
    
      IF ((R == missing) .OR. (LA + LB + LC - R) == 0 ) THEN
        ETS(k,l) = missing
      ELSE
        ETS(k,l) = REAL(LA - R)/REAL(LA + LB + LC - R)
      ENDIF
      
      IF ((LA+LC) == 0) THEN
        Bias(k,l) = missing
      ELSE
        Bias(k,l) = REAL(LA+LB)/REAL(LA+LC)
      ENDIF         
      
      CTA(k,l) = LA
      CTB(k,l) = LB
      CTC(k,l) = LC
      CTD(k,l) = LD
    END DO          
  END DO          

END SUBROUTINE calcetsb

!########################################################################
!########################################################################
!######                                                            ######
!######                  SUBROUTINE CVTPCPUNIT                     ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE cvtpcpunit(nx,ny,ntime,precip,factor,missing)
  
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Converts units of precipitation array.
!
! AUTHOR:  Eric Kemp, March 2000
!
!-----------------------------------------------------------------------
! 
! Variable declarations
! 
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nx,ny,ntime
  REAL,INTENT(INOUT) :: precip(nx,ny,ntime)
  REAL,INTENT(IN) :: factor,missing  

!-----------------------------------------------------------------------
! 
! Other variables
! 
!-----------------------------------------------------------------------

  INTEGER :: i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO k = 1,ntime
    DO j = 1,ny
      DO i = 1,nx
        IF (precip(i,j,k).ne.missing) THEN
          precip(i,j,k) = precip(i,j,k)*factor
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE cvtpcpunit
