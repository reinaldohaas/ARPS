  PROGRAM sathdf2arps
!
! Read and remap satellite data in hdf5 format from the FY2C Satellite.
! Processes the cloud cover fraction and cloud top temperature data files.
!
! Keith Brewster, CAPS
! November 1, 2007
!
! Keith Brewster, CAPS
! 04 April 2012
! Improved documentation and error handling
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: nxpix=2288
  INTEGER, PARAMETER :: nypix=2288
  REAL :: rdata(nxpix,nypix)
  REAL :: satlat(nxpix,nypix)
  REAL :: satlon(nxpix,nypix)
  REAL :: satx(nxpix,nypix)
  REAL :: saty(nxpix,nypix)
!
  CHARACTER(LEN=256) :: satinfile
  CHARACTER(LEN=256) :: sfname
  CHARACTER(LEN=12) :: satellite
  CHARACTER(LEN=6) :: satname
  CHARACTER(LEN=3) :: productid
  CHARACTER(LEN=64) :: prodname
  CHARACTER(LEN=12) :: NOMtype
  REAL :: CenterLat
  REAL :: CenterLon
  REAL :: SatHeight
  INTEGER :: CalTabYear
  INTEGER :: CalTabMonth
  INTEGER :: CalTabDay
  INTEGER :: CalTabHour
  INTEGER :: CalTabMinute
  INTEGER :: CalTabSecond
  INTEGER :: StartYear
  INTEGER :: StartMonth
  INTEGER :: StartDay
  INTEGER :: StartHour
  INTEGER :: StartMinute
  INTEGER :: StartSecond
  INTEGER :: EndYear
  INTEGER :: EndMonth
  INTEGER :: EndDay
  INTEGER :: EndHour
  INTEGER :: EndMinute
  INTEGER :: EndSecond
  INTEGER :: ProcessYear
  INTEGER :: ProcessMonth
  INTEGER :: ProcessDay
  INTEGER :: ProcessHour
  INTEGER :: ProcessMinute
  INTEGER :: ProcessSecond
  CHARACTER(LEN=128) :: comment
!
  REAL :: LwrValid,UprValid,FillValue
  REAL :: EarthRad,SampAng,StepAng,ObRec
!
! Attributes
!
  REAL :: nptsinv,pct,rscale,rbase,scpix
!
! Arps variables
!
  INTEGER :: nx,ny
  REAL :: dx,dy,xnw,ynw,ysw
  CHARACTER(LEN=80) :: arunnam
  INTEGER, ALLOCATABLE :: kntarps(:,:)
  REAL, ALLOCATABLE :: avgarps(:,:)
  REAL, ALLOCATABLE :: valarps(:,:,:)
!
! Misc. variables
!
  CHARACTER(LEN=6) :: varid
  CHARACTER(LEN=6) :: varout(2)
  CHARACTER(LEN=6) :: varunits
  CHARACTER(LEN=20) :: varname
  CHARACTER(LEN=256) :: dirname

  INTEGER :: narg
  INTEGER :: istatus
  INTEGER :: i,j,imid,jmid,iloc,jloc,hindex,di,dj
  INTEGER :: kntmiss,kntfill,i4time,kntgd
  INTEGER :: iradius,ibgn,iend,jbgn,jend,ii,jj,kntpt,isum,jsum
  INTEGER :: foutfmt,hdfcompr,mpiflag,isource,nfield
  INTEGER :: itimstr,itimend,itimsat
  INTEGER :: iyear,iyr,imon,iday,ihour,imin,isec
  INTEGER, PARAMETER :: nhist=20
  INTEGER :: knt(0:nhist)

  LOGICAL :: cldamt,tbb

  REAL, PARAMETER :: eps = 1.0E-25
  REAL, PARAMETER :: filtlen =10.0E03
  REAL :: declin,cosdec,sindec,ir0,r0,s0,s1,rlat,dlat,s1new
  REAL :: ascen,cosasc,sinasc,coslat,ealat,rl0,sl1,rlon,dlon,sl1new
  REAL :: ratio,x1,y1,x2,errval,valtmp,valmean,xpt,ypt
  REAL :: pi,dtr,rtd,dxinv,dyinv,time
  REAL :: alsq(4,4)
  REAL :: ylsq(4)
  REAL :: xlsq(4)
  REAL :: work(4,5),work1d(5)
  REAL :: sumval,avgmax,avgmin

  INCLUDE 'mp.inc'

  mp_opt=0

  pi=acos(-1.0)
  dtr=pi/180.
  rtd=180./pi

  arunnam='NULL'
  nx=103
  ny=103
  dx=3000.
  dy=3000.
  xnw=0.
  ynw=0.

  imid=(nxpix/2)+1
  jmid=(nypix/2)+1

!
! Obtain command line arguments
!
  CALL GET_COMMAND_ARGUMENT(1, satinfile, narg, istatus )
  IF( narg < 1 ) THEN
    WRITE(6,'(a)') ' Usage: sathdf52arps satellite_filename < arps.input'
    STOP
  END IF

! Read-in ARPS-grid input parameters

  CALL inisatarps(arunnam,nx,ny,dx,dy,xnw,ynw)
  WRITE(6,'(a,a)') ' Back from inisatarps runname:',TRIM(arunnam)
  WRITE(6,'(a,2i6,2f10.1)') 'nx,nx  dx,ny:', nx,ny,dx,dy
  WRITE(6,'(a,2f10.1)') ' xnw,ynw (km): ',(0.001*xnw),(0.001*ynw)
  ysw=ynw-(nx*dy)
  WRITE(6,'(a,2f10.1)') ' xnw,ysw (km): ',(0.001*xnw),(0.001*ysw)
  CALL setorig(1,xnw,ysw)
  WRITE(6,'(a)') ' lltoxy test: '
  CALL lltoxy(1,1,30.0,115.0,xpt,ypt)
  WRITE(6,'(a,2f10.1)') ' 30,115 x,y(km):',(0.001*xpt),(0.001*ypt)
!
! Allocate variables needed for remapping/averaging
!
  ALLOCATE(avgarps(nx,ny))
  ALLOCATE(valarps(nx,ny,2))
  ALLOCATE(kntarps(nx,ny))

  WRITE(6,'(2a)') 'satinfile = ',TRIM(satinfile)

  CALL readsathdf(nxpix,nypix,satinfile,rdata,                         &
     satellite,productid,prodname,NOMtype,                             &
     cldamt,tbb,CenterLat,CenterLon,SatHeight,                         &
     CalTabYear,CalTabMonth,CalTabDay,CalTabHour,                      &
     CalTabMinute,CalTabSecond,                                        &
     StartYear,StartMonth,StartDay,StartHour,StartMinute,StartSecond,  &
     EndYear,EndMonth,EndDay,EndHour,EndMinute,EndSecond,              &
     ProcessYear,ProcessMonth,ProcessDay,ProcessHour,                  &
     ProcessMinute,ProcessSecond,                                      &
     comment,LwrValid,UprValid,FillValue,EarthRad,SampAng,StepAng,     &
     ObRec,istatus)

  IF(istatus == 0) THEN

    knt=0
    kntmiss=0

    CALL ctim2abss( StartYear,StartMonth,StartDay,StartHour,           &
                    StartMinute,StartSecond,itimstr)
    CALL ctim2abss( EndYear,EndMonth,EndDay,EndHour,                   &
                    EndMinute,EndSecond,itimend)
    itimsat=itimstr+((itimend-itimstr)/2)
    CALL abss2ctim( itimsat,iyear,imon,iday,ihour,imin,isec )
    WRITE(6,'(a,i4.4,5(a,i2.2))') 'Satellite time: ',                  &
        iyear,'/',imon,'/',iday,'-',ihour,':',imin,':',isec

    IF(cldamt) THEN
      rscale=float(nhist)
      rbase=0.
      varid='satccf'
      varname='Cloud Fraction'
      varunits='None'
    ELSE
      rscale=float(nhist)/(UprValid-LwrValid)
      rbase=LwrValid
      varid='sattbb'
      varname='Brightness Temperature'
      varunits='K'
    END IF
    IF(cldamt) THEN
      scpix=255.
    ELSE
      scpix=255./(UprValid-LwrValid)
    END IF
!
!   Create histogram of valid data
!
    DO j=1,nypix
      kntgd=0
      DO i=1,nxpix
        IF( rdata(i,j) >= LwrValid .AND. rdata(i,j) <= UprValid) THEN
          hindex=INT((rdata(i,j)-rbase)*rscale)
          knt(hindex)=knt(hindex)+1
          kntgd=kntgd+1
        ELSE
          kntmiss=kntmiss+1
        END IF
      END DO
!     print *, ' Row: ',j,' Good pts: ',kntgd
    END DO

    nptsinv=100.0/(float(nxpix*nypix))
    WRITE(6,'(//a)') ' Data Histogram '
    WRITE(6,'(a)') ' Lower Limit   Count   Percent'
    DO hindex=0,20
      pct=nptsinv*knt(hindex)
      WRITE(6,'(f6.2,i12,f6.1)') (rbase+(hindex/rscale)),knt(hindex),pct
    END DO
    pct=nptsinv*kntmiss
    WRITE(6,'(/a,i12,a,f6.1,a)')                                       &
         ' Missing data count: ',kntmiss,' = ',pct,' percent'
  
    WRITE(6,'(a)') ' Calling satlatlon: '
    WRITE(6,'(a,2f12.4)') ' CenterLat, CenterLon: ',CenterLat,CenterLon
    CALL satlatlon(nxpix,nypix,CenterLat,CenterLon,                    &
                   SatHeight,EarthRad,StepAng,satlat,satlon,istatus)

    WRITE(6,'(a,2f12.4)') ' mid point lat,lon: ',                      &
          satlat(imid,jmid),satlon(imid,jmid)
!
    WRITE(6,'(a)') ' Calling lltoxy: '
    CALL lltoxy(nxpix,nypix,satlat,satlon,satx,saty)
    WRITE(6,'(a,2f12.1)') ' mid point x,y(km): ',                      &
         (0.001*satx(imid,jmid)),(0.001*saty(imid,jmid))
    WRITE(6,'(a,f12.4)') ' mid point rdata: ',rdata(imid,jmid)

    WRITE(6,'(a)') ' Beginning remap averaging/interpolation:'
    avgarps=0.
    kntarps=0
    dxinv=1.0/dx
    dyinv=1.0/dy
    DO j=1,nypix
      DO i=1,nxpix
        IF(rdata(i,j) > -1.0) THEN
          iloc=1+NINT(satx(i,j)*dxinv)
          jloc=1+NINT(saty(i,j)*dyinv)
          IF(iloc > 0 .AND. iloc <= nx .AND.                           &
             jloc > 0 .AND. jloc <= ny) THEN
            kntarps(iloc,jloc)=kntarps(iloc,jloc)+1
            avgarps(iloc,jloc)=avgarps(iloc,jloc)+rdata(i,j)
          END IF
        END IF
      END DO
    END DO
!
!  Find average from sums
!
    WRITE(6,'(a)') ' Normalization loop:'
    kntfill=0
    DO j=1,ny
      DO i=1,nx
        IF(kntarps(i,j) > 0) THEN
          avgarps(i,j)=avgarps(i,j)/float(kntarps(i,j))
        ELSE
          kntfill=kntfill+1
          avgarps(i,j)=-999.
        END IF
      END DO
    END DO
    WRITE(6,'(a,i12,a,f10.2)') ' missing points: ',kntfill,            &
          '  percent: ',(100.*float(kntfill)/float(nx*ny))
!
!   Histogram of remapped data
!
    knt=0
    kntmiss=0
    DO j=1,ny
      DO i=1,nx
        IF( avgarps(i,j) > -1.0) THEN
          hindex=INT((avgarps(i,j)-rbase)*rscale)
          hindex=max(min(hindex,20),0)
          knt(hindex)=knt(hindex)+1
        ELSE
          kntmiss=kntmiss+1
        END IF
      END DO
    END DO

    nptsinv=100.0/(float(nx*ny))
    WRITE(6,'(//a)') ' Data Histogram '
    WRITE(6,'(a)') ' Lower Limit   Count   Percent'
    DO hindex=0,20
      pct=nptsinv*knt(hindex)
      WRITE(6,'(f6.2,i12,f6.1)') (rbase+(hindex/rscale)),knt(hindex),pct
    END DO
    pct=nptsinv*kntmiss
    WRITE(6,'(/a,i12,a,f6.1,a)') &
         ' Missing data count: ',kntmiss,' = ',pct,' percent'
!
!   Use least-squares linear fit to fill any gaps.
!
    WRITE(6,'(a)') ' Least-squares gap fill:'
    DO j=1,ny
      DO i=1,nx
        IF(kntarps(i,j) > 0) THEN
          valarps(i,j,1)=avgarps(i,j)
        ELSE
          DO iradius=1,8
            kntpt=0
            sumval=0.
            alsq=0.
            ylsq=0.
            xlsq=0.
            avgmax=-999.
            avgmin=999.
            ibgn=max((i-iradius),1)
            iend=min((i+iradius),nx)
            jbgn=max((j-iradius),1)
            jend=min((j+iradius),ny)
            DO jj=jbgn,jend
              DO ii=ibgn,iend
                IF(kntarps(ii,jj) > 0) THEN
                  kntpt=kntpt+1
                  di=(ii-i)
                  dj=(jj-j)
                  avgmax=max(avgmax,avgarps(ii,jj))
                  avgmin=min(avgmin,avgarps(ii,jj))
                  sumval=sumval+avgarps(ii,jj)
                  alsq(1,1)=alsq(1,1)+di*di
                  alsq(1,2)=alsq(1,2)+dj*di
                  alsq(1,3)=alsq(1,3)+di*dj*di
                  alsq(1,4)=alsq(1,4)+di
                  alsq(2,1)=alsq(2,1)+di*dj
                  alsq(2,2)=alsq(2,2)+dj*dj
                  alsq(2,3)=alsq(2,3)+di*dj*dj
                  alsq(2,4)=alsq(2,4)+dj
                  alsq(3,1)=alsq(3,1)+di*di*dj
                  alsq(3,2)=alsq(3,2)+dj*di*dj
                  alsq(3,3)=alsq(3,3)+di*dj*di*dj
                  alsq(3,4)=alsq(3,4)+di*dj
                  alsq(4,1)=alsq(4,1)+di
                  alsq(4,2)=alsq(4,2)+dj
                  alsq(4,3)=alsq(4,3)+di*dj
                  alsq(4,4)=alsq(4,4)+1.0
                  ylsq(1)=ylsq(1)+avgarps(ii,jj)*di
                  ylsq(2)=ylsq(2)+avgarps(ii,jj)*dj
                  ylsq(3)=ylsq(3)+avgarps(ii,jj)*di*dj
                  ylsq(4)=ylsq(4)+avgarps(ii,jj)
                END IF
              END DO
            END DO
            IF(kntpt > 5) THEN
              WRITE(6,'(a,i5)') ' Found sufficient knt at iradius=',   &
                    iradius
              valmean=sumval/float(kntpt)
              CALL gjelim(4,alsq,ylsq,xlsq,work,work1d,eps,istatus)
              IF(istatus == 0) THEN
                WRITE(6,'(a,4f12.2)') ' d,valmean,avgmin,avgmax: ',    &
                       xlsq(4),valmean,avgmin,avgmax
                valarps(i,j,1)=min(max(xlsq(4),avgmin),avgmax)
                WRITE(6,'(a,f12.2)') ' Assigned valarps=',valarps(i,j,1)
              ELSE
                valarps(i,j,1)=min(max(valmean,avgmin),avgmax)
                WRITE(6,'(a,f12.2)') ' gjelim error, assigned mean=',  &
                       valarps(i,j,1)
              END IF
              EXIT
            END IF
          END DO
        END IF
      END DO
    END DO
!
! Histogram of remapped and filled data
!
    knt=0
    avgmin=999.
    avgmax=-999.
    kntmiss=0
    DO j=1,ny
      DO i=1,nx
        IF( valarps(i,j,1) > -1.0 ) THEN
          avgmin=min(avgmin,valarps(i,j,1))
          avgmax=max(avgmax,valarps(i,j,1))
          hindex=INT((valarps(i,j,1)-rbase)*rscale)
          hindex=max(min(hindex,20),0)
          knt(hindex)=knt(hindex)+1
        ELSE
          kntmiss=kntmiss+1
        END IF
      END DO
    END DO

    nptsinv=100.0/(float(nx*ny))
    WRITE(6,'(//a)') ' Data Histogram '
    WRITE(6,'(a)') ' Lower Limit   Count   Percent'
    DO hindex=0,20
      pct=nptsinv*knt(hindex)
      WRITE(6,'(f6.2,i12,f6.1)') (rbase+(hindex/rscale)),knt(hindex),pct
    END DO
    pct=nptsinv*kntmiss
    WRITE(6,'(/a,i12,a,f6.1,a)')                                         &
           ' Missing data count: ',kntmiss,' = ',pct,' percent'
    WRITE(6,'(a,f10.4,a,f10.4)') ' valmin= ',avgmin,' valmax: ',avgmax
!
! Output the array for plotting
!
    foutfmt=1
    hdfcompr=2
    mpiflag=0
    time=0.
    dirname='./'

    CALL wrtvar2(nx,ny,1,valarps,varid,varname,varunits,time,arunnam,    &
                   dirname,foutfmt,hdfcompr,mpiflag,istatus)
    IF(tbb) THEN
      nfield=2
      varout(1)=varid
      varout(2)='satcft'
      CALL coldfilt(nx,ny,filtlen,valarps(1,1,1),valarps(1,1,2))
    ELSE
      nfield=1
      varout(1)=varid
    END IF
!
! Output the array for ADAS
!
    iyr=mod(iyear,100)
    WRITE(sfname,'(a,a,i2.2,4i2.2,2a)')                                  &
       TRIM(arunnam),'.',iyr,imon,iday,ihour,imin,'.FY2C.',TRIM(varid)
    WRITE(6,'(2a)') ' sfname= ',TRIM(sfname)
    satname=satellite(1:6)
    isource=3
    CALL wtsatfld(nx,ny,nfield,                                          &
                      sfname,satname,CenterLat,CenterLon,                &
                      iyr,imon,iday,ihour,imin,isec,isource,             &
                      foutfmt,hdfcompr,                                  &
                      varout,valarps)
    WRITE(6,'(a)') ' Normal end to sathdf52arps'
  ELSE    ! bad status from readsathdf
    WRITE(6,'(a,i6)') ' Bad status from readsathdf:',istatus
    WRITE(6,'(a)') ' Abnormal end to sathdf52arps'
  END IF
  STOP
  END PROGRAM sathdf2arps

  SUBROUTINE readsathdf(nxpix,nypix,satinfile,rdata,                   &
     satellite,productid,prodname,NOMtype,                             &
     cldamt,tbb,CenterLat,CenterLon,SatHeight,                         &
     CalTabYear,CalTabMonth,CalTabDay,CalTabHour,                      &
     CalTabMinute,CalTabSecond,                                        &
     StartYear,StartMonth,StartDay,StartHour,StartMinute,StartSecond,  &
     EndYear,EndMonth,EndDay,EndHour,EndMinute,EndSecond,              &
     ProcessYear,ProcessMonth,ProcessDay,ProcessHour,                  &
     ProcessMinute,ProcessSecond,                                      &
     comment,LwrValid,UprValid,FillValue,EarthRad,SampAng,StepAng,     &
     ObRec,istatus)
!
! Keith Brewster, CAPS
! November 1, 2007
!
  USE HDF5   ! HDF5 modules
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: nxpix,nypix
  CHARACTER (LEN=256), INTENT(IN) :: satinfile
  REAL, INTENT(OUT) :: rdata
  CHARACTER(LEN=12), INTENT(OUT) :: satellite
  CHARACTER(LEN=3), INTENT(OUT) :: productid
  CHARACTER(LEN=64), INTENT(OUT):: prodname
  CHARACTER(LEN=12), INTENT(OUT):: NOMtype
  LOGICAL, INTENT(OUT):: cldamt,tbb
  REAL, INTENT(OUT):: CenterLat
  REAL, INTENT(OUT) :: CenterLon
  REAL, INTENT(OUT) :: SatHeight
  INTEGER, INTENT(OUT) :: CalTabYear
  INTEGER, INTENT(OUT) :: CalTabMonth
  INTEGER, INTENT(OUT) :: CalTabDay
  INTEGER, INTENT(OUT) :: CalTabHour
  INTEGER, INTENT(OUT) :: CalTabMinute
  INTEGER, INTENT(OUT) :: CalTabSecond
  INTEGER, INTENT(OUT) :: StartYear
  INTEGER, INTENT(OUT) :: StartMonth
  INTEGER, INTENT(OUT) :: StartDay
  INTEGER, INTENT(OUT) :: StartHour
  INTEGER, INTENT(OUT) :: StartMinute
  INTEGER, INTENT(OUT) :: StartSecond
  INTEGER, INTENT(OUT) :: EndYear
  INTEGER, INTENT(OUT) :: EndMonth
  INTEGER, INTENT(OUT) :: EndDay
  INTEGER, INTENT(OUT) :: EndHour
  INTEGER, INTENT(OUT) :: EndMinute
  INTEGER, INTENT(OUT) :: EndSecond
  INTEGER, INTENT(OUT) :: ProcessYear
  INTEGER, INTENT(OUT) :: ProcessMonth
  INTEGER, INTENT(OUT) :: ProcessDay
  INTEGER, INTENT(OUT) :: ProcessHour
  INTEGER, INTENT(OUT) :: ProcessMinute
  INTEGER, INTENT(OUT) :: ProcessSecond
  CHARACTER(LEN=128), INTENT(OUT) :: comment
  REAL, INTENT(OUT) :: LwrValid
  REAL, INTENT(OUT) :: UprValid
  REAL, INTENT(OUT) :: FillValue
  REAL, INTENT(OUT) :: EarthRad
  REAL, INTENT(OUT) :: SampAng
  REAL, INTENT(OUT) :: StepAng
  REAL, INTENT(OUT) :: ObRec
  INTEGER, INTENT(OUT) :: istatus
!
! Satellite Dataset Variables
!
  CHARACTER (LEN=256) :: LayerName
  INTEGER(HSIZE_T) :: data_dims(2)
  INTEGER(HID_T) :: file_id
  INTEGER(HID_T) :: attr1_id
  INTEGER(HID_T) :: attr2_id
  INTEGER(HID_T) :: attr3_id
  INTEGER(HID_T) :: attr4_id
  INTEGER(HID_T) :: dset_id
  INTEGER(HID_T) :: dtype_id
  INTEGER(HID_T) :: dt1_id
  INTEGER(HID_T) :: dt2_id
  INTEGER(HID_T) :: dspace_id
  INTEGER(HID_T) :: atype1_id
  INTEGER(HID_T) :: atype2_id
  INTEGER(HID_T) :: atype3_id
  INTEGER(HID_T) :: atype4_id
  INTEGER(HID_T) :: attr_type
  INTEGER(HID_T) :: aspace_id
  INTEGER(SIZE_T) :: typesize
  INTEGER(SIZE_T) :: type_size
  INTEGER(SIZE_T) :: attrlen
  INTEGER(SIZE_T) :: offset

  REAL(KIND=8) :: EA
  REAL(KIND=8) :: SamplingAngle
  REAL(KIND=8) :: SteppingAngle
  REAL(KIND=8) :: ObRecFlat
!
! Misc. Local Variables
!
  INTEGER :: iloc,lens
  CHARACTER(LEN=1) :: nullchar
!
  data_dims=1
  nullchar=ACHAR(0)

  CALL h5open_f(istatus)
!
! Open an existing file.
!
  CALL h5fopen_f (satinfile, H5F_ACC_RDONLY_F, file_id, istatus)

  IF( istatus /= 0 ) THEN
    WRITE(6,'(3a)') ' Error opening ',TRIM(satinfile),', ending...'
    RETURN
  END IF

  WRITE(6,'(a,i6)') ' Opened your file, istatus =',istatus
!
! Open dataset 
!
  cldamt=.false.
  tbb=.false.
  iloc=INDEX(satinfile,'CTA')
  IF(iloc > 0) THEN
    cldamt=.true.
    WRITE(6,'(a)') 'Looking for FY2C CLM Hourly Cloud Amount'
    CALL h5dopen_f(file_id,'FY2C CLM Hourly Cloud Amount',dset_id,istatus)
    WRITE(6,'(a,i6)') ' Opened Dataset, istatus =',istatus
    IF( istatus /= 0 ) THEN
      WRITE(6,'(a)') ' Cloud Amount Dataset not found in CTA file'
      istatus=-1
      RETURN
    END IF
  ELSE
    iloc=INDEX(satinfile,'TBB')
    IF(iloc > 0) THEN
      tbb=.true.
      WRITE(6,'(a)') 'Looking for FY2C TBB Hourly Product'
      CALL h5dopen_f(file_id,'FY2C TBB Hourly Product',dset_id,istatus)
      WRITE(6,'(a,i6)') ' Open Dataset, istatus =',istatus
      IF( istatus /= 0 ) THEN
        WRITE(6,'(a)') ' TBB Dataset not found in TBB file'
        istatus=-1
        RETURN
      END IF
    END IF
  END IF
!
! Open LayerName attribute. 
!
  CALL h5aopen_name_f (dset_id, 'LayerName', attr1_id, istatus)
  WRITE(6,'(a,i6)') ' Opened LayerName attribute, istatus =',istatus
  print *, ' attr1_id: ',attr1_id
!
! Get other attrbute space
!
  CALL h5aget_space_f(attr1_id, aspace_id, istatus)
  print *, ' aspace_id, istatus: ',aspace_id,istatus
!
! Get the string attrbute datatype
!
  CALL h5aget_type_f(attr1_id, atype1_id, istatus)
  print *, ' atype1_id, istatus: ',atype1_id,istatus
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, attr_type, istatus)
  print *, ' here 0 istatus: ',istatus
  CALL h5tset_size_f(attr_type, attrlen, istatus)
  print *, ' attrlen, istatus: ',attrlen,istatus
!
! Read LayerName attribute. 
!
  data_dims(1)=1
  CALL h5aread_f(attr1_id,atype1_id,LayerName,data_dims,istatus)
  lens=INDEX(LayerName,nullchar)
  IF(lens < 0) lens=LEN(LayerName)
  WRITE(6,'(3a,i6)') ' LayerName is =',TRIM(LayerName(1:lens)),        &
        '  istatus=',istatus
  CALL h5aclose_f(attr1_id,istatus)
!
! Find Lower Valid Range 
!
  CALL h5aopen_name_f(dset_id,'LowerValidRange',attr2_id,istatus)
  CALL h5aget_type_f(attr2_id,atype2_id,istatus)
  CALL h5aread_f(attr2_id,H5T_IEEE_F32LE,LwrValid,data_dims,istatus)
  WRITE(6,'(a,f10.2,a,i6)') ' Lower Valid Range=',LwrValid,            &
        '  istatus=',istatus
  CALL h5aclose_f(attr2_id,istatus)
!
! Find Upper Valid Range 
!
  CALL h5aopen_name_f(dset_id,'UpperValidRange',attr3_id,istatus)
  CALL h5aget_type_f(attr3_id,atype3_id,istatus)
  CALL h5aread_f(attr3_id,H5T_IEEE_F32LE,UprValid,data_dims,istatus)
  WRITE(6,'(a,f10.2,a,i6)') ' Upper Valid Range=',UprValid,            &
        '  istatus=',istatus
  CALL h5aclose_f(attr3_id,istatus)
!
! Find Lower Valid Range 
!
  CALL h5aopen_name_f(dset_id,'FillValue',attr4_id,istatus)
  CALL h5aget_type_f(attr4_id,atype4_id,istatus)
  CALL h5aread_f(attr4_id,H5T_IEEE_F32LE,FillValue,data_dims,istatus)
  WRITE(6,'(a,f12.1,a,i6)') ' Fill Value=',FillValue,                  &
        '  istatus=',istatus
  CALL h5aclose_f(attr4_id,istatus)
!
! Read-in the data
! 
  data_dims(1)=nxpix
  data_dims(2)=nypix
  CALL h5dget_type_f(dset_id, dtype_id, istatus)
  CALL h5dget_space_f(dset_id, dspace_id, istatus)
  CALL h5dread_f(dset_id, H5T_IEEE_F32LE, rdata, data_dims, istatus)
  CALL h5dclose_f(dset_id, istatus)
  CALL h5sclose_f(dspace_id, istatus)
  CALL h5tclose_f(dtype_id, istatus)
!
! Read the file info
!
  data_dims=1
  offset = 0
  CALL h5dopen_f(file_id,'NomFileInfo',dset_id,istatus)
  WRITE(6,'(a,i6)') ' Opened file info dataset, istatus:',istatus
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, dt2_id, istatus)
  typesize = 12
  CALL h5tset_size_f(dt2_id, typesize, istatus)
  CALL h5tget_size_f(dt2_id, type_size, istatus)
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "Satellite", offset, dt2_id, istatus)
  CALL h5dread_f(dset_id, dt1_id, satellite, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
  CALL h5tclose_f(dt2_id,istatus)
!
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, dt2_id, istatus)
  typesize = 3
  CALL h5tset_size_f(dt2_id, typesize, istatus)
  CALL h5tget_size_f(dt2_id, type_size, istatus)
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "ProductID", offset, dt2_id, istatus)
  CALL h5dread_f(dset_id, dt1_id, productid, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
  CALL h5tclose_f(dt2_id,istatus)

  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, dt2_id, istatus)
  typesize = 64
  CALL h5tset_size_f(dt2_id, typesize, istatus)
  CALL h5tget_size_f(dt2_id, type_size, istatus)
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "ProductName", offset, dt2_id, istatus)
  CALL h5dread_f(dset_id, dt1_id, prodname, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
  CALL h5tclose_f(dt2_id,istatus)

  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, dt2_id, istatus)
  typesize = 12
  CALL h5tset_size_f(dt2_id, typesize, istatus)
  CALL h5tget_size_f(dt2_id, type_size, istatus)
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "NOMType", offset, dt2_id, istatus)
  CALL h5dread_f(dset_id, dt1_id, NOMType, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
  CALL h5tclose_f(dt2_id,istatus)
!
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, dt2_id, istatus)
  typesize = 128
  CALL h5tset_size_f(dt2_id, typesize, istatus)
  CALL h5tget_size_f(dt2_id, type_size, istatus)
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "Comment", offset, dt2_id, istatus)
  CALL h5dread_f(dset_id, dt1_id, Comment, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
  CALL h5tclose_f(dt2_id,istatus)
!
  WRITE(6,'(2a)') 'Sat: ',TRIM(satellite)
  WRITE(6,'(2a)') 'ID:: ',TRIM(productid)
  WRITE(6,'(2a)') 'Name: ',TRIM(prodname)
  WRITE(6,'(2a)') 'NOMtype: ',TRIM(NOMtype)
  WRITE(6,'(2a)') 'Comment: ',TRIM(Comment)

  CALL h5tget_size_f(H5T_IEEE_F32LE, type_size, istatus)
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "NOMCenterLat", offset, H5T_IEEE_F32LE, istatus)
  CALL h5dread_f(dset_id, dt1_id, CenterLat, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tget_size_f(H5T_IEEE_F32LE, type_size, istatus)
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "NOMCenterLon", offset, H5T_IEEE_F32LE, istatus)
  CALL h5dread_f(dset_id, dt1_id, CenterLon, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tget_size_f(H5T_IEEE_F32LE, type_size, istatus)
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "NOMSatHeight", offset, H5T_IEEE_F32LE, istatus)
  CALL h5dread_f(dset_id, dt1_id, SatHeight, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  WRITE(6,'(a,f12.4)') 'CenterLat: ',CenterLat
  WRITE(6,'(a,f12.4)') 'CenterLon: ',CenterLon
  WRITE(6,'(a,f16.1)') 'SatHeight: ',SatHeight

  CalTabYear = 0
  CalTabMonth= 0
  CalTabDay  = 0
  CalTabHour = 0
  CalTabMinute = 0
  CalTabSecond = 0

  CALL h5tget_size_f(H5T_NATIVE_INTEGER, type_size, istatus)
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'CalTabCreateYear',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, CalTabYear, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'CalTabCreateMonth',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, CalTabMonth, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'CalTabCreateDay',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, CalTabDay, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'CalTabCreateHour',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, CalTabHour, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'CalTabCreateMinute',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, CalTabMinute, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'CalTabCreateSecond',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, CalTabSecond, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  WRITE(6,'(a,i6)') 'CalTab Year: ',CalTabYear
  WRITE(6,'(a,i6)') 'CalTab Month: ',CalTabMonth
  WRITE(6,'(a,i6)') 'CalTab Day: ',CalTabDay
  WRITE(6,'(a,i6)') 'CalTab Hour: ',CalTabHour
  WRITE(6,'(a,i6)') 'CalTab Minute: ',CalTabMinute
  WRITE(6,'(a,i6)') 'CalTab Second: ',CalTabSecond

  StartYear = 0
  StartMonth= 0
  StartDay  = 0
  StartHour = 0
  StartMinute = 0
  StartSecond = 0

  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'StartYear',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, StartYear, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'StartMonth',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, StartMonth, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'StartDay',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, StartDay, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'StartHour',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, StartHour, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'StartMinute',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, StartMinute, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'StartSecond',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, StartSecond, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  WRITE(6,'(a,i6)') 'Start Year: ',StartYear
  WRITE(6,'(a,i6)') 'Start Month: ',StartMonth
  WRITE(6,'(a,i6)') 'Start Day: ',StartDay
  WRITE(6,'(a,i6)') 'Start Hour: ',StartHour
  WRITE(6,'(a,i6)') 'Start Minute: ',StartMinute
  WRITE(6,'(a,i6)') 'Start Second: ',StartSecond

  EndYear = 0
  EndMonth= 0
  EndDay  = 0
  EndHour = 0
  EndMinute = 0
  EndSecond = 0

  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'EndYear',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, EndYear, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'EndMonth',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, EndMonth, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'EndDay',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, EndDay, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'EndHour',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, EndHour, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'EndMinute',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, EndMinute, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'EndSecond',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, EndSecond, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  WRITE(6,'(a,i6)') 'End Year: ',EndYear
  WRITE(6,'(a,i6)') 'End Month: ',EndMonth
  WRITE(6,'(a,i6)') 'End Day: ',EndDay
  WRITE(6,'(a,i6)') 'End Hour: ',EndHour
  WRITE(6,'(a,i6)') 'End Minute: ',EndMinute
  WRITE(6,'(a,i6)') 'End Second: ',EndSecond

  ProcessYear = 0
  ProcessMonth= 0
  ProcessDay  = 0
  ProcessHour = 0
  ProcessMinute = 0
  ProcessSecond = 0

  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'ProcessYear',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, ProcessYear, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'ProcessMonth',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, ProcessMonth, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'ProcessDay',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, ProcessDay, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'ProcessHour',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, ProcessHour, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'ProcessMinute',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, ProcessMinute, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F,type_size,dt1_id,istatus)
  CALL h5tinsert_f(dt1_id,'ProcessSecond',offset,H5T_NATIVE_INTEGER,istatus)
  CALL h5dread_f(dset_id, dt1_id, ProcessSecond, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  WRITE(6,'(a,i6)') 'Proc Year: ',ProcessYear
  WRITE(6,'(a,i6)') 'Proc Month: ',ProcessMonth
  WRITE(6,'(a,i6)') 'Proc Day: ',ProcessDay
  WRITE(6,'(a,i6)') 'Proc Hour: ',ProcessHour
  WRITE(6,'(a,i6)') 'Proc Minute: ',ProcessMinute
  WRITE(6,'(a,i6)') 'Proc Second: ',ProcessSecond
!
  EarthRad = 0.
  SamplingAngle = 0.
  SteppingAngle = 0.
  ObRecFlat = 0.
  CALL h5tget_size_f(H5T_NATIVE_DOUBLE, type_size, istatus)
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "EA", offset, H5T_NATIVE_DOUBLE, istatus)
  CALL h5dread_f(dset_id, dt1_id, EA, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "SamplingAngle", offset, H5T_NATIVE_DOUBLE, istatus)
  CALL h5dread_f(dset_id, dt1_id, SamplingAngle, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "SteppingAngle", offset, H5T_NATIVE_DOUBLE, istatus)
  CALL h5dread_f(dset_id, dt1_id, SteppingAngle, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dt1_id, istatus)
  CALL h5tinsert_f(dt1_id, "ObRecFlat", offset, H5T_NATIVE_DOUBLE, istatus)
  CALL h5dread_f(dset_id, dt1_id, ObRecFlat, data_dims, istatus)
  CALL h5tclose_f(dt1_id,istatus)
!
  WRITE(6,'(a,f12.1)') 'EA: ',EA
  WRITE(6,'(a,f12.1)') 'Sampling Angle: ',SamplingAngle
  WRITE(6,'(a,f12.1)') 'Stepping Angle: ',SteppingAngle
  WRITE(6,'(a,f12.1)') 'ObRecFlat: ',ObRecFlat
!
! Close the file.
!
  CALL h5fclose_f(file_id, istatus)
  WRITE(6,'(a,i6)') ' Closed your file, istatus =',istatus
!
! Close FORTRAN predefined datatypes.
  CALL h5close_f(istatus)
!
  EarthRad=EA
  StepAng=SteppingAngle
  SampAng=SamplingAngle
  ObRec=ObRecFlat
!
  RETURN
  END SUBROUTINE readsathdf
!
  SUBROUTINE satlatlon(nxpix,nypix,CenterLat,CenterLon, &
                SatHeight,EarthRad,StepAng,satlat,satlon,istatus)
!
! Calculate latitude and longitude of each satellite pixel.
!
! Keith Brewster, CAPS
! 2007
!
  IMPLICIT NONE
!
! Arguments
!
  INTEGER,INTENT(IN) :: nxpix,nypix
  REAL,INTENT(IN) :: CenterLat
  REAL,INTENT(IN) :: CenterLon
  REAL,INTENT(IN) :: SatHeight
  REAL,INTENT(IN) :: EarthRad
  REAL,INTENT(IN) :: StepAng
  REAL,INTENT(OUT) :: satlat(nxpix,nypix)
  REAL,INTENT(OUT) :: satlon(nxpix,nypix)
  INTEGER, INTENT(OUT) :: istatus
!
! Misc Local Variables
!
  INTEGER, PARAMETER :: maxiter=50
  INTEGER, PARAMETER :: istride=5
  INTEGER, PARAMETER :: jstride=5
  INTEGER :: iter,i,imid,j,jmid
  INTEGER :: i0,i1,j0,j1
  REAL :: pi,rtd
  REAL :: declin,sigi
  REAL :: c0,c1,c2,c3,c4,cx0,cx,cx2
  REAL :: a,b,c,d,z
  REAL :: clat0,clat,clat2,slat,xlat0,xlat,xlat2
  REAL :: s1,clon,xlon,xlon2
  REAL :: wistr,wjstr,wj0,wj1,wi0,wi1
  LOGICAL :: noerr
!
  pi=acos(-1.0)
  rtd=180./pi 
  wistr=1.0/float(istride)
  wjstr=1.0/float(jstride)
  c0=(SatHeight/EarthRad)+1.0
  imid=(nxpix/2)+1
  jmid=(nypix/2)+1
  satlat=-999.
  satlon=-999.
  DO j=1, nypix, jstride
    declin=(jmid-j)*StepAng
    c1=tan(declin)
    c2=c0*c1
    a=(c1*c1)+1.0
    b=-2.0*c1*c2
    c=(c2*c2)-1.0
    d=(b*b)-(4.0*a*c)
    IF(d >= 0.) THEN
      clat0=(-b+SQRT(d))/(2.0*a)
      clat0=min(max(clat0,-1.0),1.0)
      xlat0=sign(acos(clat0),declin)
!     print *, ' j: ',j,' decl= ',(rtd*declin),'  lat0=',(rtd*xlat0)
      cx0=((SatHeight/EarthRad)+1.0)/clat0
      DO i=1,nxpix,istride
        sigi=(i-imid)*StepAng
        c3=tan(sigi)
        c4=c3*cx0
        a=(c3*c3)+1.0
        b=-2.0*c3*c4
        c=(c4*c4)-1.0
        d=(b*b)-(4.0*a*c)
        IF(d >= 0.) THEN
          clon=(-b+SQRT(d))/(2.0*a)
          clon=min(max(clon,-1.0),1.0)
          xlon=sign(acos(clon),sigi)
!         print *, ' i: ',i,' sigi= ',(rtd*sigi),'  lon=',(rtd*xlon)
!
!  After obtaining an initial guess, iterate to improve solution
!
          noerr=.TRUE.
          clat=clat0
          xlat=xlat0
          DO iter=1,maxiter
            z=c1*(SatHeight+EarthRad*(1.0-clon*clat))
            slat=(z/EarthRad)
            slat=min(max(slat,-1.0),1.0)
            xlat2=sign(asin(slat),declin)
!           print *, ' iter,xlat,xlat2=',iter,(rtd*xlat),(rtd*xlat2)
            xlat=xlat2
            clat2=cos(xlat2)
            cx2=((SatHeight/EarthRad)+1.0)/clat2
            c4=c3*cx2
            b=-2.0*c3*c4
            c=(c4*c4)-1.0
            d=(b*b)-(4.0*a*c)
            IF(d >= 0.) THEN
              clon=(-b+SQRT(d))/(2.0*a)
              clon=min(max(clon,-1.0),1.0)
              xlon2=sign(acos(clon),sigi)
!             print *, ' iter,xlon,xlon2=',iter,(rtd*xlon),(rtd*xlon2)
              xlon=xlon2
            ELSE
!             print *, ' i: ',i,' sigi= ',(rtd*sigi),'  d < 0',d
              noerr=.FALSE.
              EXIT
            END IF
          END DO
          IF(noerr) THEN
            satlat(i,j)=CenterLat+(rtd*xlat)
            satlon(i,j)=CenterLon+(rtd*xlon)
          END IF
        ELSE
!         print *, ' i: ',i,' sigi= ',(rtd*sigi),'  d < 0',d
        END IF
      END DO  ! ipix
!
!  Repeat for last row
!
      IF(mod((nxpix-1),istride) /= 0) THEN
        i=nxpix
        sigi=(i-imid)*StepAng
        c3=tan(sigi)
        c4=c3*cx0
        a=(c3*c3)+1.0
        b=-2.0*c3*c4
        c=(c4*c4)-1.0
        d=(b*b)-(4.0*a*c)
        IF(d >= 0.) THEN
          clon=(-b+SQRT(d))/(2.0*a)
          clon=min(max(clon,-1.0),1.0)
          xlon=sign(acos(clon),sigi)
!         print *, ' i: ',i,' sigi= ',(rtd*sigi),'  lon=',(rtd*xlon)
!
!  After obtaining an initial guess, iterate to improve solution
!
          noerr=.TRUE.
          clat=clat0
          xlat=xlat0
          DO iter=1,maxiter
            z=c1*(SatHeight+EarthRad*(1.0-clon*clat))
            slat=(z/EarthRad)
            slat=min(max(slat,-1.0),1.0)
            xlat2=sign(asin(slat),declin)
!           print *, ' iter,xlat,xlat2=',iter,(rtd*xlat),(rtd*xlat2)
            xlat=xlat2
            clat2=cos(xlat2)
            cx2=((SatHeight/EarthRad)+1.0)/clat2
            c4=c3*cx2
            b=-2.0*c3*c4
            c=(c4*c4)-1.0
            d=(b*b)-(4.0*a*c)
            IF(d >= 0.) THEN
              clon=(-b+SQRT(d))/(2.0*a)
              clon=min(max(clon,-1.0),1.0)
              xlon2=sign(acos(clon),sigi)
!             print *, ' iter,xlon,xlon2=',iter,(rtd*xlon),(rtd*xlon2)
              xlon=xlon2
            ELSE
!             print *, ' i: ',i,' sigi= ',(rtd*sigi),'  d < 0',d
              noerr=.FALSE.
              EXIT
            END IF
          END DO
          IF(noerr) THEN
            satlat(i,j)=CenterLat+(rtd*xlat)
            satlon(i,j)=CenterLon+(rtd*xlon)
          END IF
        ELSE
!         print *, ' i: ',i,' sigi= ',(rtd*sigi),'  d < 0',d
        END IF
      END IF  ! ipix
!     print *, '  '
    ELSE
!     print *, ' decl= ',(rtd*declin),'  d < 0',d
    END IF
  END DO  ! jpix
!
! Repeat for last column
!
  IF(mod((nypix-1),jstride) /= 0) THEN
    j=nypix
    declin=(jmid-j)*StepAng
    c1=tan(declin)
    c2=c0*c1
    a=(c1*c1)+1.0
    b=-2.0*c1*c2
    c=(c2*c2)-1.0
    d=(b*b)-(4.0*a*c)
    IF(d >= 0.) THEN
      clat0=(-b+SQRT(d))/(2.0*a)
      clat0=min(max(clat0,-1.0),1.0)
      xlat0=sign(acos(clat0),declin)
!     print *, ' j: ',j,' decl= ',(rtd*declin),'  lat0=',(rtd*xlat0)
      cx0=((SatHeight/EarthRad)+1.0)/clat0
      DO i=1,nxpix,istride
        sigi=(i-imid)*StepAng
        c3=tan(sigi)
        c4=c3*cx0
        a=(c3*c3)+1.0
        b=-2.0*c3*c4
        c=(c4*c4)-1.0
        d=(b*b)-(4.0*a*c)
        IF(d >= 0.) THEN
          clon=(-b+SQRT(d))/(2.0*a)
          clon=min(max(clon,-1.0),1.0)
          xlon=sign(acos(clon),sigi)
!         print *, ' i: ',i,' sigi= ',(rtd*sigi),'  lon=',(rtd*xlon)
!
!  After obtaining an initial guess, iterate to improve solution
!
          noerr=.TRUE.
          clat=clat0
          xlat=xlat0
          DO iter=1,maxiter
            z=c1*(SatHeight+EarthRad*(1.0-clon*clat))
            slat=(z/EarthRad)
            slat=min(max(slat,-1.0),1.0)
            xlat2=sign(asin(slat),declin)
!           print *, ' iter,xlat,xlat2=',iter,(rtd*xlat),(rtd*xlat2)
            xlat=xlat2
            clat2=cos(xlat2)
            cx2=((SatHeight/EarthRad)+1.0)/clat2
            c4=c3*cx2
            b=-2.0*c3*c4
            c=(c4*c4)-1.0
            d=(b*b)-(4.0*a*c)
            IF(d >= 0.) THEN
              clon=(-b+SQRT(d))/(2.0*a)
              clon=min(max(clon,-1.0),1.0)
              xlon2=sign(acos(clon),sigi)
!             print *, ' iter,xlon,xlon2=',iter,(rtd*xlon),(rtd*xlon2)
              xlon=xlon2
            ELSE
!             print *, ' i: ',i,' sigi= ',(rtd*sigi),'  d < 0',d
              noerr=.FALSE.
              EXIT
            END IF
          END DO
          IF(noerr) THEN
            satlat(i,j)=CenterLat+(rtd*xlat)
            satlon(i,j)=CenterLon+(rtd*xlon)
          END IF
        ELSE
!         print *, ' i: ',i,' sigi= ',(rtd*sigi),'  d < 0',d
        END IF
      END DO  ! ipix
!
!  Repeat for last row
!
      IF(mod((nxpix-1),istride) /= 0) THEN
        i=nxpix
        sigi=(i-imid)*StepAng
        c3=tan(sigi)
        c4=c3*cx0
        a=(c3*c3)+1.0
        b=-2.0*c3*c4
        c=(c4*c4)-1.0
        d=(b*b)-(4.0*a*c)
        IF(d >= 0.) THEN
          clon=(-b+SQRT(d))/(2.0*a)
          clon=min(max(clon,-1.0),1.0)
          xlon=sign(acos(clon),sigi)
!         print *, ' i: ',i,' sigi= ',(rtd*sigi),'  lon=',(rtd*xlon)
!
!  After obtaining an initial guess, iterate to improve solution
!
          noerr=.TRUE.
          clat=clat0
          xlat=xlat0
          DO iter=1,maxiter
            z=c1*(SatHeight+EarthRad*(1.0-clon*clat))
            slat=(z/EarthRad)
            slat=min(max(slat,-1.0),1.0)
            xlat2=sign(asin(slat),declin)
!           print *, ' iter,xlat,xlat2=',iter,(rtd*xlat),(rtd*xlat2)
            xlat=xlat2
            clat2=cos(xlat2)
            cx2=((SatHeight/EarthRad)+1.0)/clat2
            c4=c3*cx2
            b=-2.0*c3*c4
            c=(c4*c4)-1.0
            d=(b*b)-(4.0*a*c)
            IF(d >= 0.) THEN
              clon=(-b+SQRT(d))/(2.0*a)
              clon=min(max(clon,-1.0),1.0)
              xlon2=sign(acos(clon),sigi)
!             print *, ' iter,xlon,xlon2=',iter,(rtd*xlon),(rtd*xlon2)
              xlon=xlon2
            ELSE
!             print *, ' i: ',i,' sigi= ',(rtd*sigi),'  d < 0',d
              noerr=.FALSE.
              EXIT
            END IF
          END DO
          IF(noerr) THEN
            satlat(i,j)=CenterLat+(rtd*xlat)
            satlon(i,j)=CenterLon+(rtd*xlon)
          END IF
        ELSE
!         print *, ' i: ',i,' sigi= ',(rtd*sigi),'  d < 0',d
        END IF
      END IF  ! ipix
!     print *, '  '
    ELSE
!     print *, ' decl= ',(rtd*declin),'  d < 0',d
    END IF
  END IF  ! nypix
!
! Use linear interpolation to fill-in between stride pts
!
  DO j=1,nypix
    j0=(((j-1)/jstride)*jstride)+1
    j1=j0+jstride
    IF(j1 > nypix) THEN
      j1=nypix
      wj1=float(j-j0)/float(j1-j0)
    ELSE
      wj1=wjstr*(j-j0)
    END IF
    wj0=1.0-wj1
    DO i=1,nxpix
      IF(satlat(i,j) < -90.0 ) THEN
        i0=(((i-1)/istride)*istride)+1
        i1=i0+istride
        IF(i1 > nxpix) THEN
          i1=nxpix
          wi1=float(i-i0)/float(i1-i0)
        ELSE
          wi1=wistr*(i-i0)
        END IF
        wi0=1.0-wi1
        IF(satlat(i0,j0) > -91.0 .AND. satlat(i0,j1) < 91.0 .AND.      &
           satlat(i1,j0) > -91.0 .AND. satlat(i1,j1) < 91.0 ) THEN
          satlat(i,j)=wi0*wj0*satlat(i0,j0) + wi0*wj1*satlat(i0,j1) + &
                      wi1*wj0*satlat(i1,j0) + wi1*wj1*satlat(i1,j1)
          satlon(i,j)=wi0*wj0*satlon(i0,j0) + wi0*wj1*satlon(i0,j1) + &
                      wi1*wj0*satlon(i1,j0) + wi1*wj1*satlon(i1,j1)
        END IF
      END IF
    END DO
  END DO
  istatus=0
  RETURN
  END SUBROUTINE satlatlon
! 
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE GJELIM                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  SUBROUTINE gjelim(n,a,rhs,sol,work,work1d,eps,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Solve an N x N system of equations.       [a][sol]=[rhs]
!    Using Gauss-Jordan with array pivoting.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  09/26/00
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    n        Matrix order
!    a        Array
!    rhs      Right-hand side vector of matrix equation
!    eps      Error checking threshold
!
!  OUTPUT:
!    sol      Solution vector
!    istatus  Status indicator
!
!  WORK SPACE:
!    work     Work Array
!    work1d   Work vector
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
!
  INTEGER, INTENT(IN)  :: n
  REAL,    INTENT(IN)  :: a(n,n)
  REAL,    INTENT(IN)  :: rhs(n)
  REAL,    INTENT(OUT) :: sol(n)
  REAL,    INTENT(OUT) :: work(n,n+1)
  REAL,    INTENT(OUT) :: work1d(n+1)
  REAL,    INTENT(IN)  :: eps
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!
  REAL :: pivot,const
  INTEGER np1
  INTEGER :: i,j,k,m
!
!-----------------------------------------------------------------------
!
! Initialize the work array
! First set all elements to zero.
! Fill nxn with elements from input a
! Fill last column with RHS vector.
!
!-----------------------------------------------------------------------
!
  np1=n+1

  DO j=1, np1
    DO i=1, n
      work(i,j)=0.0
    END DO
  END DO

  DO j=1, n
    DO i=1, n
      work(i,j)=a(i,j)
    END DO
  END DO

  DO i=1,n
    work(i,np1)=rhs(i)
  END DO

  DO j=1, n
!
!-----------------------------------------------------------------------
!
! Find largest element in column j
!
!-----------------------------------------------------------------------
!
    m=j
    pivot=ABS(work(m,j))
    DO i=j+1,n
      IF(ABS(work(i,j)) > pivot ) THEN
        m=i
        pivot=ABS(work(m,j))
      END IF
    END DO
!
!-----------------------------------------------------------------------
!
! Error trapping
!
!-----------------------------------------------------------------------
!
    IF( pivot < eps ) THEN
      DO i=1, n
        sol(i)=0.
      END DO
      istatus=-1
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
! Swap rows
!
!-----------------------------------------------------------------------
!
    IF(m /= j) THEN
      DO k=1,np1
        work1d(k)=work(j,k)
      END DO
      DO k=1,np1
        work(j,k)=work(m,k)
        work(m,k)=work1d(k)
      END DO
    END IF
!
!-----------------------------------------------------------------------
!
! Normalize Row
!
!-----------------------------------------------------------------------
!
    const=1./work(j,j)
    DO k=1,np1
      work(j,k)=const*work(j,k)
    END DO
    work(j,j)=1.0
!
!-----------------------------------------------------------------------
!
! Elimination
!
!-----------------------------------------------------------------------
!
    DO i=1,n
      IF ( i /= j ) THEN
        const=work(i,j)
        DO k=1,np1
          work(i,k)=work(i,k)-const*work(j,k)
        END DO
      END IF
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
! Transfer last column to sol vector
!
!-----------------------------------------------------------------------
!
  DO i=1,n
    sol(i)=work(i,n+1)
  END DO
  istatus = 1

  RETURN
  END SUBROUTINE gjelim
