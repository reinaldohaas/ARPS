!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE get_ncraddims               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
SUBROUTINE get_ncraddims(fname,indir,iradtype,ngate,nazim,nvar,          &
                         istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Obtain dimension information from a NetCDF Radar file.
!  Handles CASA Tier 2a and Tier 2b files and EEC OU-PRIME files.
!  Reports the type of file via iradtype variable.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  2007/02/27
!
!  MODIFICATION HISTORY:
!
!  07 July 2009 Keith Brewster
!  Added reading of OU-PRIME data
!
!  24 August 2010 Keith Brewster
!  Added reading of TDWR data
!
!  15 December 2010 Keith Brewster
!  Added reading of OU-PRIME Tier-2a style data
!
!  26 April 2011 Keith Brewster
!  Added reading of Foray NetCDF for field-project reprocessed data
!
!-----------------------------------------------------------------------
!
!  iradtype: 21 CASA Tier 2a Data
!            22 CASA Tier 2b Data (WDSS-II files)
!            31 OU-PRIME EEC Radar, similar to CASA Tier 2a
!            32 OU-PRIME EEC Radar, similar to CASA Tier 2b
!            42 KOUN Radar, similar to CASA Tier 2b
!            52 TDWR Radar as WDSS-II files
!            61 Foray output radar files (e.g., VORTEX2)
!            72 TDWR or NEXRAD NIDS files converted to netCDF
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=180), INTENT(IN) :: indir
  INTEGER, INTENT(OUT) :: iradtype
  INTEGER, INTENT(OUT) :: ngate
  INTEGER, INTENT(OUT) :: nazim
  INTEGER, INTENT(OUT) :: nvar
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
! Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  CHARACTER(LEN=NF_MAX_NAME) :: varname
  CHARACTER(LEN=NF_MAX_NAME) :: plugin_name
  CHARACTER(LEN=32) :: radarname
  INTEGER :: ncmode,ncid,gateid,azimid,k,lenstr
!
!-----------------------------------------------------------------------
!
!  Open netcdf file
!
!-----------------------------------------------------------------------
!
  iradtype=0
  radarname='NULL'
  cdfname = ' '
  ngate=0
  nazim=0
  nvar=0
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(fname)
  write(6,'(a,a,a)') ' Reading file dimensions from: *',TRIM(cdfname),'*'
  ncmode=nf_share
  istatus=nf_open(TRIM(cdfname),ncmode,ncid)
  IF( istatus /= NF_NOERR) THEN
    write(6,'(a,a)') ' Error opening file, message reported by netCDF = ',NF_STRERROR(istatus)
    RETURN
  END IF

  istatus=nf_inq_nvars(ncid,nvar)
  IF( istatus /= NF_NOERR) THEN
    write(6,'(a,a)') ' Error finding nvar, message reported by netCDF = ',NF_STRERROR(istatus)
    RETURN
  ELSE
    WRITE(6,'(a,i6)') ' nvar= ',nvar
  END IF
!
! Foray data refers to number of gates as maxCells
!
  istatus=nf_inq_dimid(ncid,'Gate',gateid)
  IF( istatus /= NF_NOERR) THEN
    istatus=nf_inq_dimid(ncid,'maxCells',gateid)
    IF( istatus /= NF_NOERR) THEN
      istatus=nf_inq_dimid(ncid,'gate',gateid)
      IF (istatus /= NF_NOERR) THEN
        write(6,'(a,a)') ' Error finding gate ID, message reported by netCDF = ',NF_STRERROR(istatus)
        RETURN
      ELSE
        iradtype=72
        write(6,'(a,i3)') ' File type is NEXRAD/TWDR Level-3 iradtype=',iradtype
      END IF
    ELSE
      iradtype=61
      write(6,'(a,i3)') ' File type is Foray NetCDF, iradtype=',iradtype
    END IF
  END IF

  istatus=nf_inq_dimlen(ncid,gateid,ngate)
  IF( istatus /= NF_NOERR) THEN
    write(6,'(a,a)') ' Error finding number of gates, message reported by netCDF = ',NF_STRERROR(istatus)
    RETURN
  END IF

  IF(iradtype == 61 ) THEN
    istatus=nf_inq_dimid(ncid,'Time',azimid)
    IF( istatus == NF_NOERR) THEN
      istatus=nf_inq_dimlen(ncid,azimid,nazim)
    ELSE
      write(6,'(a,a)') &
        ' Error finding number of azimuths, message reported by netCDF = ', &
          NF_STRERROR(istatus)
      RETURN
    END IF
  ELSE IF(iradtype == 72 ) THEN
    istatus=nf_inq_dimid(ncid,'azimuth',azimid)
    IF( istatus == NF_NOERR) THEN
      istatus=nf_inq_dimlen(ncid,azimid,nazim)
    ELSE
      write(6,'(a,a)') &
        ' Error finding number of azimuths, message reported by netCDF = ', &
          NF_STRERROR(istatus)
      RETURN
    END IF
  ELSE
    istatus=nf_inq_dimid(ncid,'Radial',azimid)
    IF( istatus == NF_NOERR) THEN
      iradtype=21
      write(6,'(a,i3)') ' File type is CASA Tier2a, iradtype=',iradtype
      istatus=nf_inq_dimlen(ncid,azimid,nazim)
      IF( istatus /= NF_NOERR) THEN
        write(6,'(a,a)') &
          ' Error finding number of azimuths, message reported by netCDF = ',   &
            NF_STRERROR(istatus)
        RETURN
      END IF
    ELSE
      istatus=nf_inq_dimid(ncid,'Azimuth',azimid)
      IF( istatus == NF_NOERR) THEN
        istatus=nf_inq_dimlen(ncid,azimid,nazim)
        IF( istatus /= NF_NOERR) THEN
          write(6,'(a,a)') &
            ' Error finding number of azimuths, message reported by netCDF = ',NF_STRERROR(istatus)
          RETURN
        END IF
        istatus=nf_get_att_text(ncid,NF_GLOBAL,'ConversionPlugin',plugin_name)
        IF(istatus == NF_NOERR) THEN
          IF( nvar > 6 ) THEN
            iradtype=31
            write(6,'(a,i3)') ' File type is EEC OU-PRIME 2a, iradtype=',iradtype
          ELSE
            iradtype=32
            write(6,'(a,i3)') ' File type is EEC OU-PRIME 2b, iradtype=',iradtype
          END IF
        ELSE
          istatus=nf_get_att_text(ncid,NF_GLOBAL,'radarName-value',radarname)
          IF(istatus == NF_NOERR .AND. radarname(1:4) == 'KOUN' ) THEN
            iradtype=42
            write(6,'(a,i3)') ' File type is KOUN dual-pol NEXRAD, iradtype=',iradtype
          ELSE IF(istatus == NF_NOERR .AND. radarname(1:1) /= 'K') THEN
            iradtype=52
            write(6,'(a,i3)') ' File type is WDSS-II, iradtype=',iradtype
          ELSE
            iradtype=22
            write(6,'(a,i3)') ' File type is CASA Tier2b, iradtype=',iradtype
          END IF
        END IF
      ELSE
        write(6,'(a,a)') ' Error finding azimuth ID, message reported by netCDF = ',NF_STRERROR(istatus)
        RETURN
      END IF
    END IF
  END IF

  WRITE(6,'(3(a,i6))') ' ngate= ',ngate,'  nazim= ',nazim,'  nvar= ',nvar

  istatus=NF_CLOSE(ncid)

  RETURN
END SUBROUTINE get_ncraddims
!
!########################################################################
!########################################################################
!#########                                                      #########
!#########             SUBROUTINE get_ncrad2ainfo               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
SUBROUTINE get_ncrad2ainfo(fname,indir,iradtype,nazim,nvar,tem_double,    &
                         radarname,radlat,radlon,radalt,                  &
                         itimcdf,frtime,elv,ncvarname,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Obtain information about a CASA or OU-PRIME Tier 2a NetCDF radar file,
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  2007/02/27
!
!  MODIFICATION HISTORY:
!  Keith Brewster, 12/15/2010
!  Added processing for OU-PRIME radar.  Added iradtype to argument list.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=180), INTENT(IN) :: indir
  INTEGER, INTENT(IN) :: iradtype
  INTEGER, INTENT(IN) :: nazim
  INTEGER, INTENT(IN) :: nvar
  REAL(KIND=8), INTENT(OUT) :: tem_double(nazim)
  CHARACTER(LEN=32), INTENT(OUT) :: radarname
  REAL, INTENT(OUT) :: radlat
  REAL, INTENT(OUT) :: radlon
  REAL, INTENT(OUT) :: radalt
  INTEGER, INTENT(OUT) :: itimcdf
  REAL, INTENT(OUT) :: frtime
  REAL, INTENT(OUT) :: elv
  CHARACTER(LEN=NF_MAX_NAME) :: ncvarname(nvar)
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
! Misc Internal Variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  CHARACTER (LEN=256) :: chtime
  INTEGER, PARAMETER :: itim1970=315619200
  INTEGER iyear,imonth,iday,ihour,imin,isec
  INTEGER :: ncmode,ncid,itimeid,elevid,iazim,k
  INTEGER :: latid,lonid,altid
  REAL :: elvsum,elvknt
  REAL(KIND=8) :: read_double
!
!-----------------------------------------------------------------------
!
!  Open netcdf file
!
!-----------------------------------------------------------------------
!
  frtime=0.0
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(fname)
  write(6,'(a,a)') ' Reading file information from: ',TRIM(cdfname)
  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
  IF( istatus /= NF_NOERR) THEN
    write(6,'(a,i6)') ' Error opening file, istatus= ',istatus
    RETURN
  END IF

  radarname=' '
  istatus=nf_get_att_text(ncid,nf_global,'RadarName',radarname)
  IF(istatus /= NF_NOERR) THEN
    istatus=nf_get_att_text(ncid,nf_global,'radarName-value',radarname)
    IF(istatus /= NF_NOERR) THEN
      istatus=nf_get_att_text(ncid,nf_global,'Instrument_Name',radarname)
    END IF
  END IF
  WRITE(6,'(a,a)') '  RadarName=',TRIM(radarname)
!
! For Foray data, which may be mobile radar data, we assume stationary
! and retrieve the first reported latitude, longitude and altitude.
!
  IF(iradtype == 61) THEN
    istatus=nf_inq_varid(ncid,'Latitude',latid)
    istatus=nf_get_var_double(ncid,latid,read_double)
    radlat=real(read_double)
    istatus=nf_inq_varid(ncid,'Longitude',lonid)
    istatus=nf_get_var_double(ncid,lonid,read_double)
    radlon=real(read_double)
    istatus=nf_inq_varid(ncid,'Altitude',altid)
    istatus=nf_get_var_double(ncid,altid,read_double)
    radalt=real(read_double)
  ELSE
    istatus=nf_get_att_real(ncid,nf_global,'Latitude',radlat)
    istatus=nf_get_att_real(ncid,nf_global,'Longitude',radlon)
    istatus=nf_get_att_real(ncid,nf_global,'Height',radalt)
  END IF
  WRITE(6,'(a,f10.4,a,f10.4,a,f8.1)') &
    ' Lat=',radlat,'  Lon=',radlon,'  Height=',radalt

  IF( iradtype == 31 ) THEN
    chtime=' '
    istatus=nf_get_att_text(ncid,nf_global,'Time',chtime)
    read(chtime,'(i4,5(i2))') iyear,imonth,iday,ihour,imin,isec
    CALL ctim2abss(iyear,imonth,iday,ihour,imin,isec,itimcdf)
    itimcdf=itimcdf-itim1970
  ELSE IF ( iradtype == 61 ) THEN
    istatus=nf_inq_varid(ncid,'base_time',itimeid)
    istatus=nf_get_var_int(ncid,itimeid,itimcdf)
  ELSE
    istatus=nf_inq_varid(ncid,'Time',itimeid)
    istatus=nf_get_var1_int(ncid,itimeid,1,itimcdf)
  END IF
  WRITE(6,'(a,i16,a,f10.1)') ' Time=',itimcdf,' FractionalTime=',frtime
!
! Find mean elevation that is used to sort tilts
!
  IF( iradtype == 31 ) THEN
    istatus=nf_get_att_real(ncid,nf_global,'Elevation',elv)
  ELSE IF( iradtype == 61 ) THEN
    istatus=nf_inq_varid(ncid,'Fixed_Angle',elevid)
    istatus=nf_get_var_real(ncid,elevid,elv)
  ELSE
    istatus=nf_inq_varid(ncid,'Elevation',elevid)
    istatus=nf_get_var_double(ncid,elevid,tem_double)
    elv=0.
    elvsum=0.
    elvknt=0.
    DO iazim=1,nazim
      IF(tem_double(iazim) > -5.0 .AND. tem_double(iazim) < 90.1) THEN
        elvsum=elvsum+tem_double(iazim)
        elvknt=elvknt+1.0
      END IF
    END DO
    IF(elvknt > 0.0) THEN
      elv=elvsum/elvknt
    END IF
  END IF
!
! Get variable names
!
  DO k=1,nvar
    istatus=nf_inq_varname(ncid,k,ncvarname(k))
    IF( istatus /= NF_NOERR) EXIT
  END DO
  istatus=NF_CLOSE(ncid)
  RETURN
END SUBROUTINE get_ncrad2ainfo
!
!########################################################################
!########################################################################
!#########                                                      #########
!#########             SUBROUTINE get_ncrad2binfo               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
SUBROUTINE get_ncrad2binfo(fname,indir,nazim,nvar,iradtype,              &
                         radarname,radlat,radlon,radalt,                 &
                         itimcdf,frtime,ivcp,elv,ncvarname,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Obtain information about a CASA Tier 2b, OU-PRIME NetCDF radar file,
!  or other NetCDF radar file stored with multiple variables in a
!  single file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  2007/02/26
!
!  MODIFICATION HISTORY:
!
!  26 April 2011 Keith Brewster
!  Added reading of Foray NetCDF for field-project reprocessed data
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=180), INTENT(IN) :: indir
  INTEGER, INTENT(IN) :: nazim
  INTEGER, INTENT(IN) :: nvar
  INTEGER, INTENT(IN) :: iradtype
  CHARACTER(LEN=32), INTENT(OUT) :: radarname
  REAL, INTENT(OUT) :: radlat
  REAL, INTENT(OUT) :: radlon
  REAL, INTENT(OUT) :: radalt
  INTEGER, INTENT(OUT) :: itimcdf
  REAL, INTENT(OUT) :: frtime
  INTEGER, INTENT(OUT) :: ivcp
  REAL, INTENT(OUT) :: elv
  CHARACTER(LEN=NF_MAX_NAME) :: ncvarname(nvar)
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
! Misc Internal Variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  CHARACTER (LEN=12) :: vcpstr
  INTEGER :: ncmode,ncid,itimeid,elevid,k
  REAL(KIND=8) :: dblval
  REAL :: elvsum,elvknt

  dblval=0.
  itimcdf=0
  frtime=0.
  ivcp=0
  elv=0.
  vcpstr='            '
!
!-----------------------------------------------------------------------
!
!  Open netcdf file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(fname)
  write(6,'(a,a)') ' Reading file information from: ',TRIM(cdfname)
  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
  IF( istatus /= NF_NOERR) THEN
    write(6,'(a,i6)') ' Error opening file, istatus= ',istatus
    RETURN
  END IF

  radarname=' '
  IF(iradtype == 72) THEN
    istatus=nf_get_att_text(ncid,nf_global,'ProductStation',radarname)
    radarname='T'//radarname(1:3)
    WRITE(6,'(a,a)') '  RadarName=',TRIM(radarname)
  ELSE
    istatus=nf_get_att_text(ncid,nf_global,'radarName-value',radarname)
    WRITE(6,'(a,a)') '  RadarName=',TRIM(radarname)
  END IF

  IF(iradtype == 72) THEN
    istatus=nf_get_att_real(ncid,nf_global,'RadarLatitude',radlat)
    istatus=nf_get_att_real(ncid,nf_global,'RadarLongitude',radlon)
    istatus=nf_get_att_real(ncid,nf_global,'RadarAltitude',radalt)
  ELSE
    istatus=nf_get_att_real(ncid,nf_global,'Latitude',radlat)
    istatus=nf_get_att_real(ncid,nf_global,'Longitude',radlon)
  END IF
  istatus=nf_get_att_real(ncid,nf_global,'Height',radalt)
  IF(radarname(1:3) == 'OKU' .AND. radalt < 300.) radalt=342.
  IF(radarname(1:4) == 'KOUN' .AND. radalt < 300.) radalt=380.
  WRITE(6,'(a,f10.4,a,f10.4,a,f8.1)') &
    ' Lat=',radlat,'  Lon=',radlon,'  Height=',radalt

  IF(iradtype == 72) THEN
    istatus=nf_inq_varid(ncid,'rays_time',itimeid)
    istatus=nf_get_var1_double(ncid,itimeid,1,dblval)
    itimcdf=NINT(0.001*dblval)
  ELSE
    istatus=nf_get_att_int(ncid,nf_global,'Time',itimcdf)
    istatus=nf_get_att_real(ncid,nf_global,'FractionalTime',frtime)
    WRITE(6,'(a,i16,a,f10.1)') ' Time=',itimcdf,' FractionalTime=',frtime
  END IF

  IF( iradtype == 72) THEN
    istatus=nf_get_att_int(ncid,nf_global,'VolumeCoveragePatternName',ivcp)
  ELSE
    istatus=nf_get_att_text(ncid,nf_global,'vcp-value',vcpstr)
    IF(istatus == NF_NOERR) THEN
      WRITE(6,'(a,a)') ' Read vcpstr as: ',TRIM(vcpstr)
      IF(iradtype == 32) THEN
        ivcp=1
        IF(vcpstr(1:6) == 'sector') ivcp = 2
      ELSE
        read(vcpstr,*) ivcp
      END IF
    ELSE
      WRITE(6,'(a,a)') ' No VCP String, setting to zero'
      ivcp = 0
    END IF
  END IF
  WRITE(6,'(a,i6)') ' VCP: ',ivcp

  IF( iradtype == 72) THEN
    istatus=nf_inq_varid(ncid,'elevation',elevid)
    istatus=nf_get_var1_real(ncid,elevid,1,elv)
  ELSE
    istatus=nf_get_att_real(ncid,nf_global,'Elevation',elv)
  END IF

  DO k=1,nvar
    istatus=nf_inq_varname(ncid,k,ncvarname(k))
    IF( istatus /= NF_NOERR) EXIT
  END DO
  istatus=NF_CLOSE(ncid)
  RETURN
END SUBROUTINE get_ncrad2binfo
!
!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE rdrftiltcdf                #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
!
SUBROUTINE rdrftiltcdf(nazim,ngate,fname,indir,iradtype,                &
                       refvarname,                                      &
                       rmisval,rngfval,itimcdf,frtime,initime,rfirstg,  &
                       bmwidth,azim,gtspc,refl,refqc)
!
!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads radar reflectivity data from CASA Tier 2b NetCDF
! or OU-PRIME EEC radar data files
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Keith Brewster (May 2005)
!
! MODIFICATIONS:
! 24 Mar 2007 (Keith Brewster)
! Added default value of bmwidth=1.0 if bmwidth is missing.
!
! 21 Mar 2009 (Keith Brewster)
! Added reading of integer GateFlag - refqc array.
!
! 07 July 2009 (Keith Brewster)
! Added processing of OU-PRIME EEC data.
!
! 17 July 2009 (Keith Brewster)
! Added processing of KOUN dual-pol Nexrad data.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Variable Declarations.
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nazim
  INTEGER, INTENT(IN) :: ngate
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=180), INTENT(IN) :: indir
  INTEGER, INTENT(IN) :: iradtype
  CHARACTER (LEN=80), INTENT(IN) :: refvarname
  REAL, INTENT(OUT) :: rmisval
  REAL, INTENT(OUT) :: rngfval
  INTEGER, INTENT(OUT) :: itimcdf
  REAL, INTENT(OUT) :: frtime
  INTEGER, INTENT(OUT) :: initime
  REAL, INTENT(OUT) :: rfirstg
  REAL, INTENT(OUT) :: bmwidth
!
  REAL, INTENT(OUT) :: azim(nazim)
  REAL, INTENT(OUT) :: gtspc(nazim)
  REAL, INTENT(OUT) :: refl(ngate,nazim)
  INTEGER, INTENT(OUT) :: refqc(ngate,nazim)
!
!-----------------------------------------------------------------------
!
!  netCDF variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  INTEGER :: istatus,ncid,ncmode
  INTEGER :: ipktyp,nbits
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: azmid,bmwid,elvid,gateid,gtwid,refid,rfqcid,itimeid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: lrname,ifsecs
  REAL :: rgate1,rgate2
  REAL(KIND=8) :: dblval
!
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  dblval=0.
  itimcdf=0
  frtime=0.
  rfirstg=-999.
  rmisval=-999.
  rngfval=-999.
!
!-----------------------------------------------------------------------
!
!  Open netcdf file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(fname)
  write(6,'(a,a)') ' Reading reflectivity data from: ',TRIM(cdfname)
  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
!
!-----------------------------------------------------------------------
!
!  Read Global Attributes
!
!-----------------------------------------------------------------------
!
  IF( iradtype == 72 ) THEN
    istatus=nf_inq_varid(ncid,'rays_time',itimeid)
    istatus=nf_get_var1_double(ncid,itimeid,1,dblval)
    itimcdf=NINT(0.001*dblval)
    istatus=nf_inq_varid(ncid,'gate',gateid)
    istatus=nf_get_var1_real(ncid,gateid,1,rgate1)
    istatus=nf_get_var1_real(ncid,gateid,2,rgate2)
    gtspc=rgate2-rgate1
    rfirstg=rgate1
  ELSE
    istatus=nf_get_att_int(ncid,nf_global,'Time',itimcdf)
    istatus=nf_get_att_real(ncid,nf_global,'FractionalTime',frtime)
  
    istatus=nf_get_att_real(ncid,nf_global,'RangeToFirstGate',rfirstg)
    IF(istatus /= NF_NOERR) rfirstg=200.
    istatus=nf_get_att_real(ncid,nf_global,'MissingData',rmisval)
    istatus=nf_get_att_real(ncid,nf_global,'RangeFolded',rngfval)
  END IF
  WRITE(6,'(a,i16,a,f10.1)') ' Time=',itimcdf,' FractionalTime=',frtime
  WRITE(6,'(a,f10.1)') ' Range 1st Gate=',rfirstg
  WRITE(6,'(a,f10.1,a,f10.1)') ' MissingData=',rmisval,' RangeFolded=',rngfval
  
  initime=itimcdf
  istatus=nf_get_att_real(ncid,nf_global,'InitialTime',initime)
!
!-----------------------------------------------------------------------
!
!  Read data variables
!
!-----------------------------------------------------------------------
!
  IF(iradtype == 72 ) THEN
    istatus=nf_inq_varid(ncid,'azimuth',azmid)
    istatus=nf_get_var_real(ncid,azmid,azim)
    bmwidth=1.0
  ELSE
    istatus=nf_inq_varid(ncid,'Azimuth',azmid)
    istatus=nf_get_var_real(ncid,azmid,azim)

    istatus=nf_inq_varid(ncid,'BeamWidth',bmwid)
    IF(istatus == NF_NOERR) THEN
      istatus=nf_get_var1_real(ncid,bmwid,1,bmwidth)
    ELSE
      bmwidth=1.0
    END IF

    istatus=nf_inq_varid(ncid,'GateWidth',gtwid)
    istatus=nf_get_var_real(ncid,gtwid,gtspc)
  END IF

  istatus=nf_inq_varid(ncid,TRIM(refvarname),refid)
  IF(istatus /= NF_NOERR) THEN
    istatus=nf_inq_varid(ncid,'Reflectivity',refid)
    IF(istatus /= NF_NOERR) THEN  ! EEC name is Corrected_Intensity
      istatus=nf_inq_varid(ncid,'Corrected_Intensity',refid)
      IF(istatus /= NF_NOERR) THEN  ! NIDS name is BaseReflectivity
        istatus=nf_inq_varid(ncid,'BaseReflectivity',refid)
      END IF
    END IF
  END IF

  IF(istatus == NF_NOERR) THEN
    istatus=nf_get_var_real(ncid,refid,refl)
  ELSE
    WRITE(6,'(a,i4)') &
        'Reflectivity or Corrected_Intensity not found istatus:',istatus
    RETURN
  END IF

  istatus=nf_inq_varid(ncid,'GateFlags',rfqcid)
  IF(istatus == NF_NOERR) THEN
    istatus=nf_get_var_int(ncid,rfqcid,refqc)
  ELSE
    refqc=1
  END IF
!
!-----------------------------------------------------------------------
!
! Close data file.
!
!-----------------------------------------------------------------------
!
  istatus=nf_close(ncid)

END SUBROUTINE rdrftiltcdf

!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE rdvrtiltcdf                #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
!
SUBROUTINE rdvrtiltcdf(nazim,ngate,fname,indir,iradtype,                &
                       velvarname,                                      &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       vnyquist,rfirstg,                                &
                       bmwidth,azim,gtspc,vnyq,radv,vrqc)
!
!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads radar radial velocity data from CASA Tier 2b NetCDF
! OU-PRIME EEC radar data files , KOUN NetCDF files or TDWR NetCDF files.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Keith Brewster (May 2005)
!
! MODIFICATIONS:
!
! 24 Mar 2007 (Keith Brewster)
! Added check for variable name AliasedVelocity in addition to RawVelocity
! Added default value of bmwidth=1.0 if bmwidth is missing.
!
! 21 Mar 2009 (Keith Brewster)
! Added reading of integer GateFlag - refqc array.
!
! 07 July 2009 (Keith Brewster)
! Added processing of OU-PRIME EEC data.
!
! 17 July 2009 (Keith Brewster)
! Added processing of KOUN dual-pol Nexrad data.
!
! 24 August 2010 (Keith Brewster)
! Added processing of TDWR data.
!
!------------------------------------------------------------------------
!
  IMPLICIT NONE
!
!------------------------------------------------------------------------
!
! Variable Declarations.
!
!------------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: nazim
  INTEGER, INTENT(IN) :: ngate
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=180), INTENT(IN) :: indir
  INTEGER, INTENT(IN) :: iradtype
  CHARACTER (LEN=80), INTENT(IN) :: velvarname
  REAL, INTENT(OUT) :: rmisval
  REAL, INTENT(OUT) :: rngfval
  INTEGER, INTENT(OUT) :: itimcdf
  REAL, INTENT(OUT) :: frtime
  INTEGER, INTENT(OUT) :: initime
  REAL, INTENT(OUT) :: vnyquist
  REAL, INTENT(OUT) :: rfirstg
  REAL, INTENT(OUT) :: bmwidth
!
  REAL, INTENT(OUT) :: azim(nazim)
  REAL, INTENT(OUT) :: gtspc(nazim)
  REAL, INTENT(OUT) :: vnyq(nazim)
  REAL, INTENT(OUT) :: radv(ngate,nazim)
  INTEGER, INTENT(OUT) :: vrqc(ngate,nazim)
!
!-----------------------------------------------------------------------
!
!  netCDF variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  INTEGER :: istatus,iostatus,ncid,ncmode
  INTEGER :: ipktyp,nbits
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: azmid,bmwid,elvid,gateid,gtwid,nyqid,velid,vrqcid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'

  INTEGER :: lrname,ifsecs,itimeid,rangeid
  CHARACTER(LEN=NF_MAX_NAME) vnyqtext
  REAL, PARAMETER :: vnyq_default = 38.3   ! CASA IP1
  REAL, PARAMETER :: rmischk = -300.0
  REAL, PARAMETER :: kts2ms = 0.5144444
  REAL(KIND=8) :: dblval
  REAL :: rgate1,rgate2
  INTEGER :: igate,jazim
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  itimcdf=-99
  frtime=0.
  vnyquist=-99.
  rfirstg=-999.
  rmisval=-999.
  rngfval=-999.
  initime=-99
  curtim=-999.
  vnyqtext=' '
!
!-----------------------------------------------------------------------
!
!  Open netcdf file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(fname)
  write(6,'(a,a)') ' Reading velocity data from: ',TRIM(cdfname)
  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
!
!-----------------------------------------------------------------------
!
!  Read Global Attributes
!
!-----------------------------------------------------------------------
!
  IF( iradtype == 72 ) THEN
    istatus=nf_inq_varid(ncid,'rays_time',itimeid)
    istatus=nf_get_var1_double(ncid,itimeid,1,dblval)
    itimcdf=NINT(0.001*dblval)
  ELSE
    istatus=nf_get_att_int(ncid,nf_global,'Time',itimcdf)
    istatus=nf_get_att_real(ncid,nf_global,'FractionalTime',frtime)
    WRITE(6,'(a,i16,a,f10.1)') ' Time=',itimcdf,' FractionalTime=',frtime
  END IF

  IF( iradtype == 72 ) THEN
    vnyquist=19.5
  ELSE
    istatus=nf_get_att_text(ncid,nf_global,'Nyquist_Vel-value',vnyqtext)
    IF( istatus == NF_NOERR) THEN
      READ(vnyqtext,*,iostat=iostatus) vnyquist
      IF(iostatus /= 0) vnyquist=vnyq_default
    ELSE
      istatus=nf_get_att_real(ncid,nf_global,'NyquistVelocity-value',vnyquist)
      IF(istatus /= NF_NOERR) vnyquist=vnyq_default
    END IF
  END IF
  WRITE(6,'(a,f10.2)') ' Nyquist Velocity=',vnyquist
  IF( iradtype == 72 ) THEN
    istatus=nf_inq_varid(ncid,'gate',rangeid)
    istatus=nf_get_var1_real(ncid,rangeid,1,rfirstg)
    rmisval=-32666.0
    rngfval=-999.0
  ELSE
    istatus=nf_get_att_real(ncid,nf_global,'RangeToFirstGate',rfirstg)
    IF(istatus /= NF_NOERR) rfirstg=0.
    istatus=nf_get_att_real(ncid,nf_global,'MissingData',rmisval)
    istatus=nf_get_att_real(ncid,nf_global,'RangeFolded',rngfval)
  END IF
  WRITE(6,'(a,f10.1)') ' Range 1st Gate=',rfirstg
  WRITE(6,'(a,f10.1,a,f10.1)') ' MissingData=',rmisval,' RangeFolded=',rngfval

  istatus=nf_get_att_text(ncid,nf_global,'Runname',runname)
  initime=itimcdf
  istatus=nf_get_att_int(ncid,nf_global,'InitialTime',initime)
  curtim=0.
  istatus=nf_get_att_real(ncid,nf_global,'ForecastSeconds',curtim)

  print *, ' Done reading global variables '
!
!-----------------------------------------------------------------------
!
!  Read data variables
!
!-----------------------------------------------------------------------
!
  IF( iradtype == 72 ) THEN
    istatus=nf_inq_varid(ncid,'azimuth',azmid)
    istatus=nf_get_var_real(ncid,azmid,azim)
    istatus=nf_get_att_real(ncid,nf_global,'ForecastSeconds',curtim)
    bmwidth=1.0
    istatus=nf_inq_varid(ncid,'gate',gateid)
    istatus=nf_get_var1_real(ncid,gateid,1,rgate1)
    istatus=nf_get_var1_real(ncid,gateid,2,rgate2)
    gtspc=rgate2-rgate1
    vnyq=35.0
  ELSE
    istatus=nf_inq_varid(ncid,'Azimuth',azmid)
    istatus=nf_get_var_real(ncid,azmid,azim)

    istatus=nf_inq_varid(ncid,'Beamwidth',bmwid)
    IF(istatus == NF_NOERR) THEN
      istatus=nf_get_var1_real(ncid,bmwid,1,bmwidth)
    ELSE
      istatus=nf_inq_varid(ncid,'BeamWidth',bmwid)
      IF(istatus == NF_NOERR) THEN
        istatus=nf_get_var1_real(ncid,bmwid,1,bmwidth)
      ELSE
        bmwidth=1.0
      END IF
    END IF

    istatus=nf_inq_varid(ncid,'GateWidth',gtwid)
    istatus=nf_get_var_real(ncid,gtwid,gtspc)

    istatus=nf_inq_varid(ncid,'NyquistVelocity',nyqid)
    IF(istatus == NF_NOERR) THEN
      istatus=nf_get_var_real(ncid,nyqid,vnyq)
    ELSE  ! set Nyquist velocity to the attribute Nyquist value
      vnyq=vnyquist
    END IF
  END IF

  istatus=nf_inq_varid(ncid,TRIM(velvarname),velid)
  IF(istatus /= NF_NOERR) THEN
    istatus=nf_inq_varid(ncid,'RawVelocity',velid)
    IF(istatus /= NF_NOERR) THEN
      istatus=nf_inq_varid(ncid,'AliasedVelocity',velid)
      IF(istatus /= NF_NOERR) THEN
        istatus=nf_inq_varid(ncid,'Radial_Velocity',velid)
        IF(istatus /= NF_NOERR) THEN
          istatus=nf_inq_varid(ncid,'RadialVelocity',velid)
          IF(istatus /= NF_NOERR) THEN
            istatus=nf_inq_varid(ncid,'Velocity',velid)
          END IF
        END IF
      END IF
    END IF
  END IF
  print *, ' velid=',velid
  istatus=nf_get_var_real(ncid,velid,radv)
  print *, ' istatus from get_var_real: ',istatus

  IF( iradtype == 72 ) THEN
    DO jazim=1,nazim
      DO igate=1,ngate
        IF(radv(igate,jazim) > rmischk) radv(igate,jazim)=kts2ms*radv(igate,jazim)
      END DO
    END DO
  END IF

  istatus=nf_inq_varid(ncid,'GateFlags',vrqcid)
  IF(istatus == NF_NOERR) THEN
    istatus=nf_get_var_int(ncid,vrqcid,vrqc)
  ELSE
    vrqc=1
  END IF
!
!-----------------------------------------------------------------------
!
! Close data file.
!
!-----------------------------------------------------------------------
!
  istatus=nf_close(ncid)

END SUBROUTINE rdvrtiltcdf
!
!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE rdtiltforay                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
!
SUBROUTINE rdtiltforay(nazim,ngate,fname,indir,                         &
              refvarname,velvarname,rhvvarname,                         &
              zdrvarname,kdpvarname,snrvarname,                         &
              dualpdata,rmisval,rngfval,itimcdf,frtime,vnyquist,rfirstg,&
              bmwidth,azim,gtspc,vnyq,refl,radv,rhohv,zdr,kdp,snr,      &
              tem_short2d)
!
!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads radar reflectivity and radial velocity data
! from Foray NetCDF radar data files.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Keith Brewster (April 2011)
!
! MODIFICATIONS:
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Variable Declarations.
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nazim
  INTEGER, INTENT(IN) :: ngate
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=180), INTENT(IN) :: indir
  CHARACTER (LEN=80), INTENT(IN) :: refvarname
  CHARACTER (LEN=80), INTENT(IN) :: velvarname
  CHARACTER (LEN=80), INTENT(IN) :: rhvvarname
  CHARACTER (LEN=80), INTENT(IN) :: zdrvarname
  CHARACTER (LEN=80), INTENT(IN) :: kdpvarname
  CHARACTER (LEN=80), INTENT(IN) :: snrvarname
  LOGICAL, INTENT(OUT) :: dualpdata
  REAL, INTENT(OUT) :: rmisval
  REAL, INTENT(OUT) :: rngfval
  INTEGER, INTENT(OUT) :: itimcdf
  REAL, INTENT(OUT) :: frtime
  REAL, INTENT(OUT) :: vnyquist
  REAL, INTENT(OUT) :: rfirstg
  REAL, INTENT(OUT) :: bmwidth
!
  REAL, INTENT(OUT) :: azim(nazim)
  REAL, INTENT(OUT) :: gtspc(nazim)
  REAL, INTENT(OUT) :: vnyq(nazim)
  REAL, INTENT(OUT) :: refl(ngate,nazim)
  REAL, INTENT(OUT) :: radv(ngate,nazim)
  REAL, INTENT(OUT) :: rhohv(ngate,nazim)
  REAL, INTENT(OUT) :: zdr(ngate,nazim)
  REAL, INTENT(OUT) :: kdp(ngate,nazim)
  REAL, INTENT(OUT) :: snr(ngate,nazim)
  INTEGER(KIND=2) :: tem_short2d(ngate,nazim)
!
!-----------------------------------------------------------------------
!
!  netCDF variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  INTEGER :: istatus,ncid,ncmode
  INTEGER :: ipktyp,nbits
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: itimeid,rfrstid,zhid,vrid,snrid,vnyqid,gtwid
  INTEGER :: rhvid,zdrid,kdpid
  INTEGER :: azmid,bmwid,rngid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: rmisdat = -999.0
  REAL :: sclfact,offset,gatesp
  REAL :: vrrng1,vrrngn
  REAL :: refmin,refmax
  REAL :: velmin,velmax
  REAL :: rhvmin,rhvmax
  REAL :: zdrmin,zdrmax
  REAL :: kdpmin,kdpmax
  REAL :: snrmin,snrmax
  INTEGER :: igate,jazim
  INTEGER(KIND=2) :: smisval
!
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  itimcdf=-99
  frtime=0.
  rfirstg=-999.
  rmisval=-999.
  rngfval=-999.
  sclfact=-999.
  offset=-999.
!
!-----------------------------------------------------------------------
!
!  Open netcdf file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(fname)
  write(6,'(a,a)') ' Reading refl and radv data from: ',TRIM(cdfname)
  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
!
!-----------------------------------------------------------------------
!
!  Read Data Attributes
!
!-----------------------------------------------------------------------
!
  istatus=nf_inq_varid(ncid,'base_time',itimeid)
  istatus=nf_get_var_int(ncid,itimeid,itimcdf)
  WRITE(6,'(a,i16,a,f10.1)') ' Time=',itimcdf,' FractionalTime=',frtime

  istatus=nf_inq_varid(ncid,'Range_to_First_Cell',rfrstid)
  istatus=nf_get_var_real(ncid,rfrstid,rfirstg)
  IF(istatus /= NF_NOERR) rfirstg=200.
  WRITE(6,'(a,f10.1)') ' Foray Range 1st Gate=',rfirstg
!
!-----------------------------------------------------------------------
!
!  Read data variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_inq_varid(ncid,'Azimuth',azmid)
  istatus=nf_get_var_real(ncid,azmid,azim)

  istatus=nf_inq_varid(ncid,'bm_width',bmwid)
  IF(istatus == NF_NOERR) THEN
    istatus=nf_get_var1_real(ncid,bmwid,1,bmwidth)
  ELSE
    bmwidth=1.0
  END IF
  WRITE(6,'(a,f10.1)') ' Foray Beam Width=',bmwidth

  istatus=nf_inq_varid(ncid,'Cell_Distance_Vector',rngid)
  IF( istatus == NF_NOERR ) THEN
    istatus=nf_get_var1_real(ncid,rngid,1,vrrng1)
    istatus=nf_get_var1_real(ncid,rngid,ngate,vrrngn)
    gatesp=(vrrngn-vrrng1)/float(ngate-1)
  ELSE
    istatus=nf_inq_varid(ncid,'Cell_Spacing',gtwid)
    istatus=nf_get_var_real(ncid,gtwid,gatesp)
  END IF
  WRITE(6,'(a,f10.1)') ' Gate spacing=',gatesp
  gtspc(:)=gatesp

  istatus=nf_inq_varid(ncid,'Nyquist_Velocity',vnyqid)
  istatus=nf_get_var_real(ncid,vnyqid,vnyquist)
  WRITE(6,'(a,f10.1)') ' Nyquist Velocity=',vnyquist
  vnyq(:)=vnyquist
!
!-----------------------------------------------------------------------
!
! Process Reflectivity data
!
!-----------------------------------------------------------------------
!
  refl(:,:)=rmisval
  istatus=nf_inq_varid(ncid,TRIM(refvarname),zhid)
  IF(istatus /= NF_NOERR) THEN
    istatus=nf_inq_varid(ncid,'ZH',zhid)
  END IF
  IF(istatus == NF_NOERR) THEN
    sclfact=1.0
    offset=0.0
    istatus=nf_get_att_int2(ncid,zhid,'missing_value',smisval)
    istatus=nf_get_att_real(ncid,zhid,'scale_factor',sclfact)
    istatus=nf_get_att_real(ncid,zhid,'add_offset',offset)
    WRITE(6,'(a,f10.3,a,f10.3)') ' Refl Scale Factor=',sclfact, &
            '  Offset=',offset

    istatus=nf_get_var_int2(ncid,zhid,tem_short2d)

    DO jazim=1,nazim
      refmin=999.
      refmax=-999.
      DO igate=1,ngate
        IF(tem_short2d(igate,jazim) /= smisval) THEN
          refl(igate,jazim)=offset+(sclfact*tem_short2d(igate,jazim))
          IF(refl(igate,jazim) > -90.) THEN
            refmin=min(refmin,refl(igate,jazim))
            refmax=max(refmax,refl(igate,jazim))
          END IF
        END IF
      END DO
!     WRITE(6,'(a,f7.1,a,f7.1,a,f7.1)') ' Azim:',azim(jazim), &
!           ' Refmin:',refmin,'  Refmax:',refmax
    END DO
  ELSE
    WRITE(6,'(a,i4)') &
      'Reflectivity not found istatus:',istatus
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
! Process Radial Velocity data
!
!-----------------------------------------------------------------------
!
  radv(:,:)=rmisval
  istatus=nf_inq_varid(ncid,TRIM(velvarname),vrid)
  IF(istatus /= NF_NOERR) THEN
    istatus=nf_inq_varid(ncid,'VE',vrid)
  END IF
  IF(istatus == NF_NOERR) THEN
    sclfact=1.0
    offset=0.0
    istatus=nf_get_att_int2(ncid,vrid,'missing_value',smisval)
    istatus=nf_get_att_real(ncid,vrid,'scale_factor',sclfact)
    istatus=nf_get_att_real(ncid,vrid,'add_offset',offset)
    WRITE(6,'(a,f10.3,a,f10.3)') ' Velocity Scale Factor=',sclfact, &
          '  Offset=',offset

    istatus=nf_get_var_int2(ncid,vrid,tem_short2d)

    DO jazim=1,nazim
      velmin=999.
      velmax=-999.
      DO igate=1,ngate
        IF(tem_short2d(igate,jazim) /= smisval) THEN
          radv(igate,jazim)=offset+(sclfact*tem_short2d(igate,jazim))
          velmin=min(velmin,radv(igate,jazim))
          velmax=max(velmax,radv(igate,jazim))
        END IF
      END DO
!     WRITE(6,'(a,f7.1,a,f7.1,a,f7.1)') ' Azim:',azim(jazim), &
!         ' Velmin:',velmin,'  Velmax:',velmax
    END DO
  ELSE
    WRITE(6,'(a,i4)') &
    'Radial Velocity (VE) not found istatus:',istatus
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
! Process RhoHV, if available
! Here set rhohv to 1.0 as default rather than a missing value.
!
!-----------------------------------------------------------------------
!
  rhohv(:,:)=1.0
  istatus=nf_inq_varid(ncid,TRIM(rhvvarname),rhvid)
  IF(istatus /= NF_NOERR) THEN
    istatus=nf_inq_varid(ncid,'RhoHV',rhvid)
  END IF
  IF(istatus == NF_NOERR) THEN
    sclfact=1.0
    offset=0.0
    istatus=nf_get_att_int2(ncid,rhvid,'missing_value',smisval)
    istatus=nf_get_att_real(ncid,rhvid,'scale_factor',sclfact)
    istatus=nf_get_att_real(ncid,rhvid,'add_offset',offset)
    WRITE(6,'(a,f10.3,a,f10.3)') ' RhoHV Scale Factor=',sclfact, &
          '  Offset=',offset

    istatus=nf_get_var_int2(ncid,rhvid,tem_short2d)

    DO jazim=1,nazim
      rhvmin=999.
      rhvmax=-999.
      DO igate=1,ngate
        IF(tem_short2d(igate,jazim) /= smisval) THEN
          rhohv(igate,jazim)=offset+(sclfact*tem_short2d(igate,jazim))
          rhvmin=min(rhvmin,rhohv(igate,jazim))
          rhvmax=max(rhvmax,rhohv(igate,jazim))
        END IF
      END DO
      WRITE(6,'(a,f7.1,a,f7.1,a,f7.1)') ' Azim:',azim(jazim), &
          ' RhoHVmin:',rhvmin,'  RhoHVmax:',rhvmax
    END DO
  ELSE
    WRITE(6,'(a,i4)') &
    'RhoHV (RhoHV) not found istatus:',istatus
    WRITE(6,'(a,f10.1)') 'Using RhoHV default:',rhohv(1,1)
  END IF
!
!-----------------------------------------------------------------------
!
! Process Zdr, if available
!
!-----------------------------------------------------------------------
!
  zdr(:,:)=rmisdat
  istatus=nf_inq_varid(ncid,TRIM(zdrvarname),zdrid)
  IF(istatus /= NF_NOERR) THEN
    istatus=nf_inq_varid(ncid,'ZDR',zdrid)
  END IF
  IF(istatus == NF_NOERR) THEN
    dualpdata=.TRUE.
    sclfact=1.0
    offset=0.0
    istatus=nf_get_att_int2(ncid,zdrid,'missing_value',smisval)
    istatus=nf_get_att_real(ncid,zdrid,'scale_factor',sclfact)
    istatus=nf_get_att_real(ncid,zdrid,'add_offset',offset)
    WRITE(6,'(a,f10.3,a,f10.3)') ' Zdr Scale Factor=',sclfact, &
          '  Offset=',offset

    istatus=nf_get_var_int2(ncid,zdrid,tem_short2d)

    DO jazim=1,nazim
      zdrmin=999.
      zdrmax=-999.
      DO igate=1,ngate
        IF(tem_short2d(igate,jazim) /= smisval) THEN
          zdr(igate,jazim)=offset+(sclfact*tem_short2d(igate,jazim))
          zdrmin=min(zdrmin,zdr(igate,jazim))
          zdrmax=max(zdrmax,zdr(igate,jazim))
        END IF
      END DO
      WRITE(6,'(a,f7.1,a,f7.1,a,f7.1)') ' Azim:',azim(jazim), &
          ' Zdrmin:',zdrmin,'  Zdrmax:',zdrmax
    END DO
  ELSE
    WRITE(6,'(a,i4)') &
    'Zdr not found istatus:',istatus
    WRITE(6,'(a,f10.1)') 'Setting Zdr to missing:',zdr(1,1)
  END IF
!
!-----------------------------------------------------------------------
!
! Process Kdp, if available
!
!-----------------------------------------------------------------------
!
  kdp(:,:)=rmisdat
  istatus=nf_inq_varid(ncid,TRIM(kdpvarname),kdpid)
  IF(istatus /= NF_NOERR) THEN
    istatus=nf_inq_varid(ncid,'Kdp',kdpid)
  END IF
  IF(istatus == NF_NOERR) THEN
    dualpdata=.TRUE.
    sclfact=1.0
    offset=0.0
    istatus=nf_get_att_int2(ncid,kdpid,'missing_value',smisval)
    istatus=nf_get_att_real(ncid,kdpid,'scale_factor',sclfact)
    istatus=nf_get_att_real(ncid,kdpid,'add_offset',offset)
    WRITE(6,'(a,f10.3,a,f10.3)') ' Kdp Scale Factor=',sclfact, &
          '  Offset=',offset

    istatus=nf_get_var_int2(ncid,kdpid,tem_short2d)

    DO jazim=1,nazim
      kdpmin=999.
      kdpmax=-999.
      DO igate=1,ngate
        IF(tem_short2d(igate,jazim) /= smisval) THEN
          kdp(igate,jazim)=offset+(sclfact*tem_short2d(igate,jazim))
          kdpmin=min(kdpmin,kdp(igate,jazim))
          kdpmax=max(kdpmax,kdp(igate,jazim))
        END IF
      END DO
      WRITE(6,'(a,f7.1,a,f7.1,a,f7.1)') ' Azim:',azim(jazim), &
          ' Kdpmin:',kdpmin,'  Kdpmax:',kdpmax
    END DO
  ELSE
    WRITE(6,'(a,i4)') &
    'Kdp (KDP) not found istatus:',istatus
    WRITE(6,'(a,f10.1)') 'Using Kdp default:',kdp(1,1)
  END IF
!
!-----------------------------------------------------------------------
!
! Process Signal-to-Noise
!
!-----------------------------------------------------------------------
!
  snr(:,:)=10.0
  istatus=nf_inq_varid(ncid,TRIM(snrvarname),snrid)
  IF( istatus /= NF_NOERR ) THEN
    istatus=nf_inq_varid(ncid,'SN',snrid)
  END IF
  IF(istatus == NF_NOERR) THEN
    sclfact=1.0
    offset=0.0
    istatus=nf_get_att_int2(ncid,snrid,'missing_value',smisval)
    istatus=nf_get_att_real(ncid,snrid,'scale_factor',sclfact)
    istatus=nf_get_att_real(ncid,snrid,'add_offset',offset)
    WRITE(6,'(a,f10.3,a,f10.3)') ' SNR Scale Factor=',sclfact, &
          '  Offset=',offset

    istatus=nf_get_var_int2(ncid,snrid,tem_short2d)

    DO jazim=1,nazim
      snrmin=999.
      snrmax=-999.
      DO igate=1,ngate
        IF(tem_short2d(igate,jazim) /= smisval) THEN
          snr(igate,jazim)=offset+(sclfact*tem_short2d(igate,jazim))
          snrmin=min(snrmin,snr(igate,jazim))
          snrmax=max(snrmax,snr(igate,jazim))
        END IF
      END DO
!     WRITE(6,'(a,f7.1,a,f7.1,a,f7.1)') ' Azim:',azim(jazim), &
!         ' SNRmin:',snrmin,'  SNRmax:',rhvmax
    END DO
  ELSE
    WRITE(6,'(a,i4)') 'SNR (SN) not found istatus:',istatus
    WRITE(6,'(a,f10.1)') 'Using SNR default:',snr(1,1)
  END IF
!
!-----------------------------------------------------------------------
!
! Close data file.
!
!-----------------------------------------------------------------------
!
  istatus=nf_close(ncid)
  RETURN

END SUBROUTINE rdtiltforay
!
!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE rdpftiltcdf                #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
!
SUBROUTINE rdpftiltcdf(nazim,ngate,fname,indir,idp,                     &
                       zvvvarname,zdrvarname,zdpvarname,                &
                       rhvvarname,kdpvarname,                           &
                       rmisval,rngfval,itimcdf,frtime,initime,rfirstg,  &
                       bmwidth,azim,gtspc,refl,refqc)
!
!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads polarimetric radar data from KOUN radar data files
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Keith Brewster (May 2005)
!
! MODIFICATIONS:
! 24 Mar 2007 (Keith Brewster)
! Added default value of beamw=1.0 if beamw is missing.
!
! 21 Mar 2009 (Keith Brewster)
! Added reading of integer GateFlag - refqc array.
!
! 07 July 2009 (Keith Brewster)
! Added processing of OU-PRIME EEC data.
!
! 17 July 2009 (Keith Brewster)
! Added processing of KOUN dual-pol Nexrad data.
!
! 01 Sept 2009 (Youngsun Jung)
! Modified for KOUN polarimetric variables.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Variable Declarations.
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nazim
  INTEGER, INTENT(IN) :: ngate
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=80), INTENT(IN) :: indir
  INTEGER, INTENT(IN) :: idp                 ! idex for dualpol variables

  CHARACTER(LEN=80), INTENT(IN) :: zvvvarname
  CHARACTER(LEN=80), INTENT(IN) :: zdrvarname
  CHARACTER(LEN=80), INTENT(IN) :: zdpvarname
  CHARACTER(LEN=80), INTENT(IN) :: rhvvarname
  CHARACTER(LEN=80), INTENT(IN) :: kdpvarname

  REAL, INTENT(OUT) :: rmisval
  REAL, INTENT(OUT) :: rngfval
  INTEGER, INTENT(OUT) :: itimcdf
  REAL, INTENT(OUT) :: frtime
  INTEGER, INTENT(OUT) :: initime
  REAL, INTENT(OUT) :: rfirstg
!
  REAL, INTENT(OUT) :: azim(nazim)
  REAL, INTENT(OUT) :: bmwidth
  REAL, INTENT(OUT) :: gtspc(nazim)
  REAL, INTENT(OUT) :: refl(ngate,nazim)
  INTEGER, INTENT(OUT) :: refqc(ngate,nazim)
!
!-----------------------------------------------------------------------
!
!  netCDF variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  INTEGER :: istatus,ncid,ncmode
  INTEGER :: ipktyp,nbits
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: nazimid,ngateid
  INTEGER :: azmid,bmwid,elvid,gtwid,refid,rfqcid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: lrname,ifsecs
!
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  itimcdf=-99
  frtime=-999.
  rfirstg=-999.
  rmisval=-999.
  rngfval=-999.
!
!-----------------------------------------------------------------------
!
!  Open netcdf file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(fname)
  write(6,'(a,a)') ' Reading polarimetric data from: ',TRIM(cdfname)
  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
!
!-----------------------------------------------------------------------
!
!  Read Global Attributes
!
!-----------------------------------------------------------------------
!
  istatus=nf_get_att_int(ncid,nf_global,'Time',itimcdf)
  istatus=nf_get_att_real(ncid,nf_global,'FractionalTime',frtime)
  WRITE(6,'(a,i16,a,f10.1)') ' Time=',itimcdf,' FractionalTime=',frtime

  istatus=nf_get_att_real(ncid,nf_global,'RangeToFirstGate',rfirstg)
  IF(istatus /= NF_NOERR) rfirstg=200.
  istatus=nf_get_att_real(ncid,nf_global,'MissingData',rmisval)
  istatus=nf_get_att_real(ncid,nf_global,'RangeFolded',rngfval)
  WRITE(6,'(a,f10.1)') ' Range 1st Gate=',rfirstg
  WRITE(6,'(a,f10.1,a,f10.1)') ' MissingData=',rmisval,' RangeFolded=',rngfval

  initime=itimcdf
  istatus=nf_get_att_real(ncid,nf_global,'InitialTime',initime)
!
!-----------------------------------------------------------------------
!
!  Read data variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_inq_varid(ncid,'Azimuth',azmid)
  istatus=nf_get_var_real(ncid,azmid,azim)

  istatus=nf_inq_varid(ncid,'BeamWidth',bmwid)
  IF(istatus == NF_NOERR) THEN
    istatus=nf_get_var1_real(ncid,bmwid,1,bmwidth)
  ELSE
    bmwidth=1.0
  END IF

  istatus=nf_inq_varid(ncid,'GateWidth',gtwid)
  istatus=nf_get_var_real(ncid,gtwid,gtspc)

  SELECT CASE (idp)
  CASE(3)
    istatus=nf_inq_varid(ncid,TRIM(zvvvarname),refid)
    IF( istatus /= NF_NOERR ) &
      istatus=nf_inq_varid(ncid,'Zvv',refid)
  CASE(4)
    istatus=nf_inq_varid(ncid,TRIM(zdrvarname),refid)
    IF( istatus /= NF_NOERR ) &
      istatus=nf_inq_varid(ncid,'Zdr',refid)
  CASE(5)
    istatus=nf_inq_varid(ncid,TRIM(zdpvarname),refid)
    IF( istatus /= NF_NOERR ) &
      istatus=nf_inq_varid(ncid,'Zdp',refid)
  CASE(6)
    istatus=nf_inq_varid(ncid,TRIM(kdpvarname),refid)
    IF( istatus /= NF_NOERR ) &
      istatus=nf_inq_varid(ncid,'Kdp',refid)
  CASE(7)
    istatus=nf_inq_varid(ncid,TRIM(rhvvarname),refid)
    IF( istatus /= NF_NOERR ) &
      istatus=nf_inq_varid(ncid,'RhoHV',refid)
  END SELECT

  IF(istatus == NF_NOERR) THEN
    istatus=nf_get_var_real(ncid,refid,refl)
  END IF

  istatus=nf_inq_varid(ncid,'GateFlags',rfqcid)
  IF(istatus == NF_NOERR) THEN
    istatus=nf_get_var_int(ncid,rfqcid,refqc)
  ELSE   ! set qc flags to "good"
    refqc=1
  END IF
!
!-----------------------------------------------------------------------
!
! Close data file.
!
!-----------------------------------------------------------------------
!
  istatus=nf_close(ncid)

END SUBROUTINE rdpftiltcdf
!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE rdvvtiltcdf                #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE rdvvtiltcdf(nazim,ngate,fname,indir,                         &
                       rmisval,rngfval,itimcdf,frtime,initime,rfirstg,  &
                       bmwidth,azim,gtspc,vort)
!
!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads vertical vorticity data from NetCDF data files
! Interpolated vertical vorticity is a special dataset generated by the
! radar emulator to use for verification of radar algortihms.
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Keith Brewster (May 2005)
!
! MODIFICATIONS:
! 24 Mar 2007 (Keith Brewster)
! Added default value of bmwidth=1.0 if bmwidth is missing.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Variable Declarations.
!
!------------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: nazim
  INTEGER, INTENT(IN) :: ngate
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=180), INTENT(IN) :: indir
  REAL, INTENT(OUT) :: rmisval
  REAL, INTENT(OUT) :: rngfval
  INTEGER, INTENT(OUT) :: itimcdf
  REAL, INTENT(OUT) :: frtime
  INTEGER, INTENT(OUT) :: initime
  REAL, INTENT(OUT) :: rfirstg
  REAL, INTENT(OUT) :: bmwidth
!
  REAL, INTENT(OUT) :: azim(nazim)
  REAL, INTENT(OUT) :: gtspc(nazim)
  REAL, INTENT(OUT) :: vort(ngate,nazim)
!
!-----------------------------------------------------------------------
!
!  netCDF variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  INTEGER :: istatus,ncid,ncmode
  INTEGER :: ipktyp,nbits
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: azmid,bmwid,elvid,gtwid,vorid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: lrname,ifsecs
!
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  itimcdf=-99
  frtime=-999.
  rfirstg=-999.
  rmisval=-999.
  rngfval=-999.
!
!-----------------------------------------------------------------------
!
!  Open netcdf file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(fname)
  write(6,'(a,a)') ' Reading vorticity data from: ',TRIM(cdfname)
  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
!
!-----------------------------------------------------------------------
!
!  Read Global Attributes
!
!-----------------------------------------------------------------------
!
  istatus=nf_get_att_int(ncid,nf_global,'Time',itimcdf)
  istatus=nf_get_att_real(ncid,nf_global,'FractionalTime',frtime)
  istatus=nf_get_att_real(ncid,nf_global,'RangeToFirstGate',rfirstg)
  IF(istatus /= NF_NOERR) rfirstg=200.
  istatus=nf_get_att_real(ncid,nf_global,'MissingData',rmisval)
  istatus=nf_get_att_real(ncid,nf_global,'RangeFolded',rngfval)

  istatus=nf_get_att_text(ncid,nf_global,'Runname',runname)
  initime=itimcdf
  istatus=nf_get_att_int(ncid,nf_global,'InitialTime',initime)
  curtim=0.
  istatus=nf_get_att_real(ncid,nf_global,'ForecastSeconds',curtim)
!
!-----------------------------------------------------------------------
!
!  Read data variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_inq_varid(ncid,'Azimuth',azmid)
  istatus=nf_get_var_real(ncid,azmid,azim)
  istatus=nf_inq_varid(ncid,'BeamWidth',bmwid)
  IF(istatus == NF_NOERR) THEN
    istatus=nf_get_var1_real(ncid,bmwid,1,bmwidth)
  ELSE
    bmwidth=1.0
  END IF
  istatus=nf_inq_varid(ncid,'GateWidth',gtwid)
  istatus=nf_get_var_real(ncid,gtwid,gtspc)
  istatus=nf_inq_varid(ncid,'VerticalVorticity',vorid)
  istatus=nf_get_var_real(ncid,vorid,vort)
!
!-----------------------------------------------------------------------
!
! Close data file.
!
!-----------------------------------------------------------------------
!
  istatus=nf_close(ncid)

END SUBROUTINE rdvvtiltcdf
!
!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE rdrhotiltcdf               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
!
SUBROUTINE rdrhotiltcdf(nazim,ngate,fname,indir,rhvvarname,             &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       rfirstg,azim,rhohv)
!
!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads radar dual-pol Rho-HV data from CASA Level-2b NetCDF
! or OU-PRIME EEC radar data files
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Keith Brewster (Nov 2006)
!
! MODIFICATIONS:
!
! 07 July 2009 (Keith Brewster)
! Added processing of OU-PRIME EEC data.
!
! 17 July 2009 (Keith Brewster)
! Added processing of KOUN dual-pol Nexrad data.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Variable Declarations.
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nazim
  INTEGER, INTENT(IN) :: ngate
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=180), INTENT(IN) :: indir
  CHARACTER (LEN=80), INTENT(IN) :: rhvvarname
  REAL, INTENT(OUT) :: rmisval
  REAL, INTENT(OUT) :: rngfval
  INTEGER, INTENT(OUT) :: itimcdf
  REAL, INTENT(OUT) :: frtime
  INTEGER, INTENT(OUT) :: initime
  REAL, INTENT(OUT) :: rfirstg
!
  REAL, INTENT(OUT) :: azim(nazim)
  REAL, INTENT(OUT) :: rhohv(ngate,nazim)
!
!-----------------------------------------------------------------------
!
!  netCDF variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  INTEGER :: istatus,ncid,ncmode
  INTEGER :: ipktyp,nbits
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: azmid,elvid,rhvid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: lrname,ifsecs
!
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  itimcdf=-99
  frtime=-999.
  rfirstg=-999.
  rmisval=-999.
  rngfval=-999.
!
!-----------------------------------------------------------------------
!
!  Open netcdf file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(fname)
  write(6,'(a,/a)') ' Reading Rho-HV data from: ',TRIM(cdfname)
  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
!
!-----------------------------------------------------------------------
!
!  Read Global Attributes
!
!-----------------------------------------------------------------------
!
  istatus=nf_get_att_int(ncid,nf_global,'Time',itimcdf)
  istatus=nf_get_att_real(ncid,nf_global,'FractionalTime',frtime)
  WRITE(6,'(a,i16,a,f10.1)') ' Time=',itimcdf,' FractionalTime=',frtime

  istatus=nf_get_att_real(ncid,nf_global,'RangeToFirstGate',rfirstg)
  IF(istatus /= NF_NOERR) rfirstg=200.
  istatus=nf_get_att_real(ncid,nf_global,'MissingData',rmisval)
  istatus=nf_get_att_real(ncid,nf_global,'RangeFolded',rngfval)
  WRITE(6,'(a,f10.1)') ' Range 1st Gate=',rfirstg
  WRITE(6,'(a,f10.1,a,f10.1)') ' MissingData=',rmisval,' RangeFolded=',rngfval

  initime=itimcdf
  istatus=nf_get_att_real(ncid,nf_global,'InitialTime',initime)
!
!-----------------------------------------------------------------------
!
!  Read data variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_inq_varid(ncid,'Azimuth',azmid)
  istatus=nf_get_var_real(ncid,azmid,azim)

  istatus=nf_inq_varid(ncid,'RhoHV',rhvid)
  IF(istatus /= NF_NOERR ) THEN
    istatus=nf_inq_varid(ncid,'RhoHV',rhvid)
  END IF
  IF(istatus == NF_NOERR ) THEN
    istatus=nf_get_var_real(ncid,rhvid,rhohv)
  ELSE
    WRITE(6,'(a,i6)') ' Error finding RhoHV Variable, istatus=',istatus
  END IF
!
!-----------------------------------------------------------------------
!
! Close data file.
!
!-----------------------------------------------------------------------
!
  istatus=nf_close(ncid)

END SUBROUTINE rdrhotiltcdf
!
!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE rd2atiltcdf                #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########    Center for Analysis and Prediction of Storms      #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
!
SUBROUTINE rd2atiltcdf(nazim,ngate,iradtype,fname,indir,                & 
                       refvarname,velvarname,rhvvarname,snrvarname,     &
                       zdrvarname,kdpvarname,                           &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       vnyquist,rfirstg,bmwidth,istrgate,               &
                       gcfst,azim,gtspc,attref,refl,radv,rhohv,         &
                       zdr,kdp,gateqc,tem_double,tem_real)
!
!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads radar data from CASA Level-2a NetCDF data files
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Keith Brewster, CAPS  (Feb 2007)
!
! MODIFICATIONS:
!
! 12/12/2011 (Y. Wang)
! Fixed bug when both global attribute "NyquistVelocity" and variable
! "Nyquist_Velocity" are missed from the CASA data file.
!
! 01/15/2013 (K. Brewster)
! Updated for CASA data file changes.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Variable Declarations.
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nazim
  INTEGER, INTENT(IN) :: ngate
  INTEGER, INTENT(IN) :: iradtype
  CHARACTER (LEN=256), INTENT(IN) :: fname
  CHARACTER (LEN=180), INTENT(IN) :: indir
  CHARACTER (LEN=80), INTENT(IN) :: refvarname
  CHARACTER (LEN=80), INTENT(IN) :: velvarname
  CHARACTER (LEN=80), INTENT(IN) :: rhvvarname
  CHARACTER (LEN=80), INTENT(IN) :: snrvarname
  CHARACTER (LEN=80), INTENT(IN) :: zdrvarname
  CHARACTER (LEN=80), INTENT(IN) :: kdpvarname
  REAL, INTENT(OUT) :: rmisval
  REAL, INTENT(OUT) :: rngfval
  INTEGER, INTENT(OUT) :: itimcdf
  REAL, INTENT(OUT) :: frtime
  INTEGER, INTENT(OUT) :: initime
  REAL, INTENT(OUT) :: vnyquist
  REAL, INTENT(OUT) :: rfirstg
  REAL, INTENT(OUT) :: bmwidth
  INTEGER, INTENT(OUT) :: istrgate(nazim)
  INTEGER, INTENT(OUT) :: gcfst(nazim)
!
  REAL, INTENT(OUT) :: azim(nazim)
  REAL, INTENT(OUT) :: gtspc(nazim)
  REAL, INTENT(OUT) :: attref(ngate,nazim)
  REAL, INTENT(OUT) :: refl(ngate,nazim)
  REAL, INTENT(OUT) :: radv(ngate,nazim)
  REAL, INTENT(OUT) :: rhohv(ngate,nazim)
  REAL, INTENT(OUT) :: zdr(ngate,nazim)
  REAL, INTENT(OUT) :: kdp(ngate,nazim)

  INTEGER, INTENT(OUT) :: gateqc(ngate,nazim)
  REAL(KIND=8), INTENT(OUT) :: tem_double(nazim)
  REAL, INTENT(OUT) :: tem_real(ngate,nazim)
!
!-----------------------------------------------------------------------
!
!  netCDF variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  CHARACTER (LEN=256) :: chtime
  INTEGER :: istatus,ncid,ncmode
  INTEGER :: ipktyp,nbits
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: itimeid,istrid,strid,vnyqid
  INTEGER :: zdrid,kdpid
  INTEGER :: bmwid,gcfid,azmid,gtwid,arfid,refid,velid,rhvid,gtqcid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: mm2m = 0.001
  REAL, PARAMETER :: vnyq_default = 38.3   ! CASA IP1
  INTEGER, PARAMETER :: itim1970=315619200
  REAL :: elvsum,rfrsum,strange
  INTEGER :: i,j
  INTEGER :: lrname,ifsecs
  INTEGER :: iyear,imonth,iday,ihour,imin,isec
  INTEGER :: jazim,kntelv,kntrfr
  INTEGER :: irfirstg,igtwid
!
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  itimcdf=-99
  frtime=-999.
  vnyquist=-99.
  rfirstg=-999.
  istrgate=1
  gcfst=1
  rmisval=-999.
  rngfval=-999.
  bmwidth=-99.
  gateqc=0

  refid=-99
  velid=-99
  rhvid=-99
  zdrid=-99
  kdpid=-99
!
!-----------------------------------------------------------------------
!
!  Open netcdf file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(indir),'/',TRIM(fname)
  write(6,'(a,/a)') ' Reading Tier 2a ref, vel, rhohv data from: ',TRIM(cdfname)
  ncmode=nf_share
  istatus=nf_open(cdfname,ncmode,ncid)
!
!-----------------------------------------------------------------------
!
!  Read Attributes, Often Global
!
!-----------------------------------------------------------------------
!
  IF( iradtype == 31 ) THEN
    chtime=' '
    istatus=nf_get_att_text(ncid,nf_global,'Time',chtime)
    read(chtime,'(i4,5(i2))') iyear,imonth,iday,ihour,imin,isec
    CALL ctim2abss(iyear,imonth,iday,ihour,imin,isec,itimcdf)
    itimcdf=itimcdf-itim1970
  ELSE
    istatus=nf_inq_varid(ncid,'Time',itimeid)
    istatus=nf_get_var1_int(ncid,itimeid,1,itimcdf)
  END IF
  frtime=0.
  initime=itimcdf
  WRITE(6,'(a,i16,a,f10.2)') ' Time=',itimcdf,'  FrTime=',frtime

  vnyquist=vnyq_default
  IF ( iradtype == 31 ) THEN
    istatus=nf_get_att_real(ncid,nf_global,'NyquistVelocity-value',vnyquist)
    istatus=nf_inq_varid(ncid,'Beamwidth',bmwid)
    istatus=nf_get_var1_real(ncid,bmwid,1,bmwidth)
  ELSE
    istatus=nf_get_att_real(ncid,nf_global,'NyquistVelocity',vnyquist)
    IF(istatus /= NF_NOERR) THEN
      istatus=nf_inq_varid(ncid,'Nyquist_Velocity',vnyqid)
      IF(istatus == NF_NOERR) THEN
        istatus=nf_get_var_real(ncid,vnyqid,vnyquist)
      ELSE
        WRITE(6,'(1x,2a,/,10x,a,f5.2,/,10x,a)') 'WARNING: Global attribute ',  &
        '"NyquistVelocity" and variable "Nyquist_Velocity" not found.',&
        'vnyquist is temporarily set to default value:',vnyquist, &
        'Set fixed_nyqv and related variables in the input file.'
      END IF
    END IF
    istatus=nf_get_att_real(ncid,nf_global,'AntennaBeamwidth',bmwidth)
  END IF
  WRITE(6,'(1x,a,f5.2)') 'Nyquist Velocity: ',vnyquist
  WRITE(6,'(1x,a,f10.2)') '  Beamwidth=',bmwidth

  IF( iradtype == 31 ) THEN
    istrgate=1
    istatus=nf_get_att_real(ncid,nf_global,'RangeToFirstGate',rfirstg)
  ELSE
    istatus=nf_inq_varid(ncid,'StartRange',strid)
    IF( istatus == NF_NOERR) THEN
      istatus=nf_get_var_double(ncid,strid,tem_double)
      rfrsum=0.
      kntrfr=0
      DO jazim=1,nazim
        strange=mm2m*tem_double(jazim)
!       print *, ' jazim: ',jazim,' strange(m): ',strange
        IF(abs(strange) < 10000.0) THEN
          rfrsum=rfrsum+strange
          kntrfr=kntrfr+1
        END IF
      END DO
      rfirstg=0.
      IF(kntrfr > 0) THEN
        rfirstg=rfrsum/float(kntrfr)
      END IF
    ELSE
      istatus=nf_inq_varid(ncid,'Range_to_First_Cell',istrid)
      print *, ' Range to 1st var: ',istrid,istatus
      IF( istatus == NF_NOERR) THEN
        istatus=nf_get_var_real(ncid,istrid,rfirstg)
        print *, ' Range to 1st gate: ',rfirstg,istatus
      ELSE
        print *, ' Getting range to first gate from variable attribute'
        istatus=nf_inq_varid(ncid,TRIM(refvarname),istrid)
        IF( istatus == NF_NOERR) THEN
          istatus=nf_get_att_int(ncid,istrid,'meters_to_first_cell',irfirstg)
          IF( istatus == NF_NOERR ) rfirstg = float(irfirstg)
          print *, ' Range to 1st gate: ',rfirstg,istatus
        END IF
      END IF
    END IF
  END IF
  WRITE(6,'(a,f10.1)') 'range to 1st gate (m)=',rfirstg

  rmisval=99900.
  rngfval=99901.
  istatus=nf_get_att_real(ncid,nf_global,'MissingData',rmisval)
  istatus=nf_get_att_real(ncid,nf_global,'RangeFolded',rngfval)
  print *, '  MissingData=',rmisval,' RangeFolded=',rngfval

!
!-----------------------------------------------------------------------
!
!  Read data variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_inq_varid(ncid,'GcfState',gcfid)
  IF(istatus == nf_noerr) THEN
    WRITE(6,'(a,i6)') ' gcfid=',gcfid
    istatus=nf_get_var_int(ncid,gcfid,gcfst)
  ELSE
    WRITE(6,'(a)') ' No Ground Clutter Filter State, assuming filter on'
  END IF

  WRITE(6,'(a,i6)') ' rd2atiltcdf: iradtype =',iradtype
  IF( iradtype == 31 ) THEN
    istatus=nf_inq_varid(ncid,'Azimuth',azmid)
    print *, ' azmid:',azmid
    istatus=nf_get_var_real(ncid,azmid,azim)

    istatus=nf_inq_varid(ncid,'GateWidth',gtwid)
    print *, ' gtwid:',gtwid
    istatus=nf_get_var_real(ncid,gtwid,gtspc)

    istatus=nf_inq_varid(ncid,TRIM(refvarname),refid)
    IF( istatus /= NF_NOERR ) THEN
      istatus=nf_inq_varid(ncid,'Corrected_Intensity',refid)
    END IF
    print *, ' refid:',refid
    istatus=nf_get_var_real(ncid,refid,refl)
    attref=refl

    istatus=nf_inq_varid(ncid,TRIM(rhvvarname),rhvid)
    IF( istatus == NF_NOERR) THEN
      WRITE(6,'(1x,a,i5)') 'Found Rho-HV, id=',rhvid
      istatus=nf_get_var_real(ncid,rhvid,rhohv)
    ELSE
      istatus=nf_inq_varid(ncid,'RhoHV',rhvid)
      IF( istatus == NF_NOERR) THEN
        WRITE(6,'(1x,a,i5)') 'Found Rho-HV, id=',rhvid
        istatus=nf_get_var_real(ncid,rhvid,rhohv)
      ELSE
        WRITE(6,'(1x,a)') 'No Rho-HV data in file'
      END IF
    END IF

    istatus=nf_inq_varid(ncid,TRIM(velvarname),velid)
    IF( istatus /= NF_NOERR ) THEN
      istatus=nf_inq_varid(ncid,'Radial_Velocity',velid)
    END IF
    print *, ' velid:',velid
    istatus=nf_get_var_real(ncid,velid,radv)

    istatus=nf_inq_varid(ncid,'Second_Trip_Echo_Mask',gtqcid)
    print *, ' gtqcid:',gtqcid
    IF(istatus == NF_NOERR) THEN
      istatus=nf_get_var_real(ncid,gtqcid,tem_real)
      DO j=1,nazim
        DO i=1,ngate
          IF( tem_real(i,j) == 0.0) gateqc(i,j)=1
        END DO
      END DO
    ELSE   ! set qc flags to "good"
      gateqc=1
    END IF

    istatus=nf_inq_varid(ncid,TRIM(zdrvarname),zdrid)
    IF( istatus == NF_NOERR) THEN
      WRITE(6,'(1x,a,i5)') 'Found Zdr, id=',zdrid
      istatus=nf_get_var_real(ncid,zdrid,zdr)
    ELSE
      istatus=nf_inq_varid(ncid,'DifferentialReflectivity',zdrid)
      IF( istatus == NF_NOERR) THEN
        WRITE(6,'(1x,a,i5)') 'Found Zdr, id=',zdrid
        istatus=nf_get_var_real(ncid,zdrid,zdr)
      ELSE
        WRITE(6,'(1x,a)') 'No Zdr data in file'
      END IF
    END IF

    istatus=nf_inq_varid(ncid,TRIM(kdpvarname),zdrid)
    IF( istatus == NF_NOERR) THEN
      WRITE(6,'(1x,a,i5)') 'Found Kdp, id=',kdpid
      istatus=nf_get_var_real(ncid,kdpid,kdp)
    ELSE
      istatus=nf_inq_varid(ncid,'Kdp',kdpid)
      IF( istatus == NF_NOERR) THEN
        WRITE(6,'(1x,a,i5)') 'Found Kdp, id=',kdpid
        istatus=nf_get_var_real(ncid,kdpid,kdp)
      ELSE
        WRITE(6,'(1x,a)') 'No Kdp data in file'
      END IF
    END IF

  ELSE

    istatus=nf_inq_varid(ncid,'Azimuth',azmid)
    print *, ' azmid:',azmid
    istatus=nf_get_var_double(ncid,azmid,tem_double)
    azim=-999.
    DO jazim=1,nazim
      azim(jazim)=tem_double(jazim)
    END DO

    istatus=nf_inq_varid(ncid,'GateWidth',gtwid)
    IF( istatus == NF_NOERR ) THEN
      print *, ' gtwid:',gtwid,istatus
      istatus=nf_get_var_double(ncid,gtwid,tem_double)
      gtspc(1)=mm2m*tem_double(1)
      DO jazim=2,nazim
        gtspc(jazim)=gtspc(1)
      END DO
    ELSE
      istatus=nf_inq_varid(ncid,'Cell Spacing',gtwid)
      print *, ' gtwid,istatus:',gtwid,istatus
      IF( istatus == NF_NOERR ) THEN
        istatus=nf_get_var_real(ncid,gtwid,gtspc(1))
        DO jazim=2,nazim
          gtspc(jazim)=gtspc(1)
        END DO
      ELSE
        print *, ' Getting gate spacing from variable attribute'
        istatus=nf_inq_varid(ncid,TRIM(refvarname),gtwid)
        IF( istatus == NF_NOERR) THEN
          istatus=nf_get_att_int(ncid,gtwid,'meters_between_cells',igtwid)
          print *, ' meters between cells: ',igtwid
          IF( istatus == NF_NOERR ) THEN
            gtspc(1)=float(igtwid)
            DO jazim=2,nazim
              gtspc(jazim)=gtspc(1)
            END DO
          END IF
        END IF
      END IF
    END IF

    istatus=nf_inq_varid(ncid,'Reflectivity',arfid)
    print *, ' arfid:',arfid
    istatus=nf_get_var_real(ncid,arfid,attref)

    istatus=nf_inq_varid(ncid,TRIM(refvarname),refid)
    IF( istatus /= NF_NOERR ) THEN
      istatus=nf_inq_varid(ncid,'CorrectedReflectivity',refid)
    END IF
    print *, ' refid:',refid
    istatus=nf_get_var_real(ncid,refid,refl)

    istatus=nf_inq_varid(ncid,TRIM(velvarname),velid)
    IF( istatus /= NF_NOERR ) THEN
      istatus=nf_inq_varid(ncid,'Velocity',velid)
    END IF
    print *, ' velid:',velid
    istatus=nf_get_var_real(ncid,velid,radv)

    istatus=nf_inq_varid(ncid,TRIM(rhvvarname),rhvid)
    IF( istatus == NF_NOERR) THEN
      WRITE(6,'(1x,a,i5)') 'Found Rho-HV, id=',rhvid
      istatus=nf_get_var_real(ncid,rhvid,rhohv)
    ELSE
      istatus=nf_inq_varid(ncid,'CrossPolCorrelation',rhvid)
      IF( istatus == NF_NOERR) THEN
        WRITE(6,'(1x,a,i5)') 'Found Rho-HV, id=',rhvid
        istatus=nf_get_var_real(ncid,rhvid,rhohv)
      ELSE
        WRITE(6,'(1x,a)') 'No Rho-HV data in file'
      END IF
    END IF

    istatus=nf_inq_varid(ncid,TRIM(zdrvarname),zdrid)
    IF( istatus == NF_NOERR) THEN
      WRITE(6,'(1x,a,i5)') 'Found Zdr, id=',zdrid
      istatus=nf_get_var_real(ncid,zdrid,zdr)
    ELSE
      istatus=nf_inq_varid(ncid,'DifferentialReflectivity',zdrid)
      IF( istatus == NF_NOERR) THEN
        WRITE(6,'(1x,a,i5)') 'Found Zdr, id=',zdrid
        istatus=nf_get_var_real(ncid,zdrid,zdr)
      ELSE
        WRITE(6,'(1x,a)') 'No Zdr data in file'
      END IF
    END IF

    istatus=nf_inq_varid(ncid,TRIM(kdpvarname),kdpid)
    IF( istatus == NF_NOERR) THEN
      WRITE(6,'(1x,a,i5)') 'Found Kdp, id=',kdpid
      istatus=nf_get_var_real(ncid,kdpid,kdp)
    ELSE
      istatus=nf_inq_varid(ncid,'Kdp',kdpid)
      IF( istatus == NF_NOERR) THEN
        WRITE(6,'(1x,a,i5)') 'Found Kdp, id=',kdpid
        istatus=nf_get_var_real(ncid,kdpid,kdp)
      ELSE
        WRITE(6,'(1x,a)') 'No Kdp data in file'
      END IF
    END IF

    istatus=nf_inq_varid(ncid,'GateFlags',gtqcid)
    IF(istatus == NF_NOERR) THEN
      istatus=nf_get_var_int(ncid,gtqcid,gateqc)
    ELSE   ! set qc flags to "good"
      gateqc=1
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
! Close data file.
!
!-----------------------------------------------------------------------
!
  istatus=nf_close(ncid)
  RETURN

END SUBROUTINE rd2atiltcdf
