SUBROUTINE wtrftiltcdf(ngate,maxazim,nazim,fname,outdir,varname,       &
                       radname,radlat,radlon,radelv,ivcp,elv,          &
                       rmisval,rngfval,itimcdf,frtime,initime,         &
                       vnyquist,rfrgate,                               &
                       azim,beamw,gtspc,refl)
!
! PURPOSE: Write one tilt of radar reflectivity data in NetCDF format.
! AUTHOR: Keith Brewster, CAPS
!
  IMPLICIT NONE
  INTEGER :: ngate
  INTEGER :: maxazim
  INTEGER :: nazim
  CHARACTER (LEN=256) :: fname
  CHARACTER (LEN=80)  :: outdir
  CHARACTER (LEN=40 ) :: varname
  CHARACTER (LEN=4)   :: radname
  REAL :: radlat
  REAL :: radlon
  REAL :: radelv
  INTEGER :: ivcp
  REAL :: elv
  REAL :: rmisval
  REAL :: rngfval
  INTEGER :: itimcdf
  REAL :: frtime
  INTEGER :: initime
  REAL :: vnyquist
  REAL :: rfrgate
!
  REAL :: azim(maxazim)
  REAL :: beamw(maxazim)
  REAL :: gtspc(maxazim)
  REAL :: refl(ngate,maxazim)
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
!  netCDF output variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: swpdim(2)
  INTEGER :: azmdim(1)
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: nazimid,ngateid
  INTEGER :: azmid,bmwid,elvid,gtwid,refid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: istr1d(1),iknt1d(1)
  INTEGER :: istr2d(2),iknt2d(2)
  INTEGER :: lrname,lvarnam,ifsecs
!
  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ncmode=nf_clobber
  lrname=LEN_TRIM(runname)
  lvarnam=LEN_TRIM(varname)
  istr1d(1)=1
  istr2d(1)=1
  istr2d(2)=1
  iknt1d(1)=nazim
  iknt2d(1)=ngate
  iknt2d(2)=nazim
!
!-----------------------------------------------------------------------
!
!  Open netCDF file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(outdir),'/',TRIM(fname)
  write(6,'(a,a)') ' Writing data to: ',TRIM(cdfname)
  istatus=nf_create(cdfname,ncmode,ncid)
!
!-----------------------------------------------------------------------
!
!  Define Dimensions
!
!-----------------------------------------------------------------------
!
  istatus=nf_def_dim(ncid,'Azimuth',nazim,nazimid)
  istatus=nf_def_dim(ncid,'Gate',ngate,ngateid)
  azmdim(1)=nazimid
  swpdim(1)=ngateid
  swpdim(2)=nazimid
!
!-----------------------------------------------------------------------
!
!  Define Variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_def_var(ncid,'Azimuth',nf_real,1,azmdim,azmid)
  istatus=nf_def_var(ncid,'BeamWidth',nf_real,1,azmdim,bmwid)
  istatus=nf_def_var(ncid,'GateWidth',nf_real,1,azmdim,gtwid)
  istatus=nf_def_var(ncid,TRIM(varname),nf_real,2,swpdim,refid)
!
!-----------------------------------------------------------------------
!
!  Set Variable Attributes
!
!-----------------------------------------------------------------------
!
  istatus=nf_put_att_text(ncid,azmid,'Units',7,'Degrees')
  istatus=nf_put_att_text(ncid,bmwid,'Units',7,'Degrees')
  istatus=nf_put_att_text(ncid,gtwid,'Units',6,'Meters')
  istatus=nf_put_att_text(ncid,refid,'Units',3,'dBZ')
!
!-----------------------------------------------------------------------
!
!  Set Global Attributes
!
!-----------------------------------------------------------------------
!
  istatus=nf_put_att_text(ncid,nf_global,'TypeName',lvarnam,TRIM(varname))
  istatus=nf_put_att_text(ncid,nf_global,'DataType',9,'RadialSet')
  istatus=nf_put_att_real(ncid,nf_global,'Latitude',nf_real,1,radlat)
  istatus=nf_put_att_real(ncid,nf_global,'Longitude',nf_real,1,radlon)
  istatus=nf_put_att_real(ncid,nf_global,'Height',nf_real,1,radelv)
  istatus=nf_put_att_int(ncid,nf_global,'Time',nf_int,1,itimcdf)
  istatus=nf_put_att_real(ncid,nf_global,'FractionalTime',nf_real,1,frtime)
  istatus=nf_put_att_text(ncid,nf_global,'attributes',26,  &
                            ' Nyquist_Vel radarName vcp')
  istatus=nf_put_att_text(ncid,nf_global,'Nyquist_Vel-unit',15,'MetersPerSecond')
  istatus=nf_put_att_real(ncid,nf_global,'Nyquist_Vel-value',nf_real,1,vnyquist)
  istatus=nf_put_att_text(ncid,nf_global,'radarName-unit',13,'dimensionless')
  istatus=nf_put_att_text(ncid,nf_global,'radarName-value',4,radname)
  istatus=nf_put_att_text(ncid,nf_global,'vcp-unit',13,'dimensionless')
  istatus=nf_put_att_int(ncid,nf_global,'vcp-value',nf_int,1,ivcp)
  istatus=nf_put_att_real(ncid,nf_global,'Elevation',nf_real,1,elv)
  istatus=nf_put_att_text(ncid,nf_global,'ElevationUnits',7,'Degrees')
  istatus=nf_put_att_real(ncid,nf_global,'RangeToFirstGate',nf_real,1,rfrgate)
  istatus=nf_put_att_text(ncid,nf_global,'RangeToFirstGateUnits',6,'Meters')
  istatus=nf_put_att_real(ncid,nf_global,'MissingData',nf_real,1,rmisval)
  istatus=nf_put_att_real(ncid,nf_global,'RangeFolded',nf_real,1,rngfval)

  istatus=nf_put_att_text(ncid,nf_global,'Runname',lrname,runname)
  istatus=nf_put_att_int(ncid,nf_global,'InitialTime',nf_int,1,initime)
  istatus=nf_put_att_real(ncid,nf_global,'ForecastSeconds',nf_real,1,curtim)
  
  istatus=nf_enddef(ncid)
!
!-----------------------------------------------------------------------
!
!  Write data variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_put_vara_real(ncid,azmid,istr1d,iknt1d,azim)
  istatus=nf_put_vara_real(ncid,bmwid,istr1d,iknt1d,beamw)
  istatus=nf_put_vara_real(ncid,gtwid,istr1d,iknt1d,gtspc)
  istatus=nf_put_vara_real(ncid,refid,istr2d,iknt2d,refl)
!
!-----------------------------------------------------------------------
!
!  Close and compress reflectivity file
!
!-----------------------------------------------------------------------
!
  istatus=nf_close(ncid)
  CALL cmprsgz(cdfname)
END SUBROUTINE wtrftiltcdf

SUBROUTINE wtvrtiltcdf(ngate,maxazim,nazim,fname,outdir,varname,       &
                       radname,radlat,radlon,radelv,ivcp,elv,          &
                       rmisval,rngfval,itimcdf,frtime,initime,         &
                       vnyquist,rfrgate,                               &
                       azim,beamw,gtspc,vnyq,radv)
!
! PURPOSE: Write one tilt of radar reflectivity data in NetCDF format.
! AUTHOR: Keith Brewster, CAPS
!
  IMPLICIT NONE
  INTEGER :: ngate
  INTEGER :: maxazim
  INTEGER :: nazim
  CHARACTER (LEN=256) :: fname
  CHARACTER (LEN=80)  :: outdir
  CHARACTER (LEN=40)  :: varname
  CHARACTER (LEN=4)   :: radname
  REAL :: radlat
  REAL :: radlon
  REAL :: radelv
  INTEGER :: ivcp
  REAL :: elv
  REAL :: rmisval
  REAL :: rngfval
  INTEGER :: itimcdf
  REAL :: frtime
  INTEGER :: initime
  REAL :: vnyquist
  REAL :: rfrgate

  REAL :: azim(maxazim)
  REAL :: beamw(maxazim)
  REAL :: gtspc(maxazim)
  REAL :: vnyq(maxazim)
  REAL :: radv(ngate,maxazim)
!
!-----------------------------------------------------------------------
!
!  netCDF variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: cdfname
  INTEGER :: istatus,ncid,ncmode
!
!-----------------------------------------------------------------------
!
!  netCDF output variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: swpdim(2)
  INTEGER :: azmdim(1)
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: nazimid,ngateid
  INTEGER :: azmid,bmwid,elvid,gtwid,nyqid,refid,velid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: itim1970=315619200
  INTEGER :: istr1d(1),iknt1d(1)
  INTEGER :: istr2d(2),iknt2d(2)
  INTEGER :: lrname,lvarnam,ifsecs,initimeS
!
  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ncmode=nf_clobber
  lrname=LEN_TRIM(runname)
  lvarnam=LEN_TRIM(varname)
  istr1d(1)=1
  istr2d(1)=1
  istr2d(2)=1
  iknt1d(1)=nazim
  iknt2d(1)=ngate
  iknt2d(2)=nazim
!
!-----------------------------------------------------------------------
!
!  Open netCDF file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(outdir),'/',TRIM(fname)
  write(6,'(a,a)') ' Writing data to: ',TRIM(cdfname)
  istatus=nf_create(cdfname,ncmode,ncid)
!
!-----------------------------------------------------------------------
!
!  Define Dimensions
!
!-----------------------------------------------------------------------
!
  istatus=nf_def_dim(ncid,'Azimuth',nazim,nazimid)
  istatus=nf_def_dim(ncid,'Gate',ngate,ngateid)
  azmdim(1)=nazimid
  swpdim(1)=ngateid
  swpdim(2)=nazimid
!
!-----------------------------------------------------------------------
!
!  Define Variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_def_var(ncid,'Azimuth',nf_real,1,azmdim,azmid)
  istatus=nf_def_var(ncid,'BeamWidth',nf_real,1,azmdim,bmwid)
  istatus=nf_def_var(ncid,'GateWidth',nf_real,1,azmdim,gtwid)
  istatus=nf_def_var(ncid,'NyquistVelocity',nf_real,1,azmdim,nyqid)
  istatus=nf_def_var(ncid,TRIM(varname),nf_real,2,swpdim,velid)
!
!-----------------------------------------------------------------------
!
!  Set Variable Attributes
!
!-----------------------------------------------------------------------
!
  istatus=nf_put_att_text(ncid,azmid,'Units',7,'Degrees')
  istatus=nf_put_att_text(ncid,bmwid,'Units',7,'Degrees')
  istatus=nf_put_att_text(ncid,gtwid,'Units',6,'Meters')
  istatus=nf_put_att_text(ncid,nyqid,'Units',15,'MetersPerSecond')
  istatus=nf_put_att_text(ncid,velid,'Units',15,'MetersPerSecond')
!
!-----------------------------------------------------------------------
!
!  Set Global Attributes
!
!-----------------------------------------------------------------------
!
  istatus=nf_put_att_text(ncid,nf_global,'TypeName',lvarnam,TRIM(varname))
  istatus=nf_put_att_text(ncid,nf_global,'DataType',9,'RadialSet')
  istatus=nf_put_att_real(ncid,nf_global,'Latitude',nf_real,1,radlat)
  istatus=nf_put_att_real(ncid,nf_global,'Longitude',nf_real,1,radlon)
  istatus=nf_put_att_real(ncid,nf_global,'Height',nf_real,1,radelv)
  istatus=nf_put_att_int(ncid,nf_global,'Time',nf_int,1,itimcdf)
  istatus=nf_put_att_real(ncid,nf_global,'FractionalTime',nf_real,1,frtime)
  istatus=nf_put_att_text(ncid,nf_global,'attributes',26,  &
                          ' Nyquist_Vel radarName vcp')
  istatus=nf_put_att_text(ncid,nf_global,'Nyquist_Vel-unit',15,'MetersPerSecond')
  istatus=nf_put_att_real(ncid,nf_global,'Nyquist_Vel-value',nf_real,1,vnyquist)
  istatus=nf_put_att_text(ncid,nf_global,'radarName-unit',13,'dimensionless')
  istatus=nf_put_att_text(ncid,nf_global,'radarName-value',4,radname)
  istatus=nf_put_att_text(ncid,nf_global,'vcp-unit',13,'dimensionless')
  istatus=nf_put_att_int(ncid,nf_global,'vcp-value',nf_int,1,ivcp)
  istatus=nf_put_att_real(ncid,nf_global,'Elevation',nf_real,1,elv)
  istatus=nf_put_att_text(ncid,nf_global,'ElevationUnits',7,'Degrees')
  istatus=nf_put_att_real(ncid,nf_global,'RangeToFirstGate',nf_real,1,rfrgate)
  istatus=nf_put_att_text(ncid,nf_global,'RangeToFirstGateUnits',6,'Meters')
  istatus=nf_put_att_real(ncid,nf_global,'MissingData',nf_real,1,rmisval)
  istatus=nf_put_att_real(ncid,nf_global,'RangeFolded',nf_real,1,rngfval)

  istatus=nf_put_att_text(ncid,nf_global,'Runname',lrname,runname)
  istatus=nf_put_att_int(ncid,nf_global,'InitialTime',nf_int,1,initime)
  istatus=nf_put_att_real(ncid,nf_global,'ForecastSeconds',nf_real,1,curtim)

  istatus=nf_enddef(ncid)
!
!-----------------------------------------------------------------------
!
!  Write data variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_put_vara_real(ncid,azmid,istr1d,iknt1d,azim)
  istatus=nf_put_vara_real(ncid,bmwid,istr1d,iknt1d,beamw)
  !istatus=nf_put_vara_real(ncid,gtwid,istr1d,iknt1d,gtspc)
  !istatus=nf_put_vara_real(ncid,gtwid,istr1d,iknt1d,vnyq)
  istatus=nf_put_vara_real(ncid,nyqid,istr1d,iknt1d,gtspc) ! Bug fix from Alex Schenkman
  istatus=nf_put_vara_real(ncid,velid,istr1d,iknt1d,vnyq)  ! on 01/19/2012
  istatus=nf_put_vara_real(ncid,gtwid,istr2d,iknt2d,radv)
!
!-----------------------------------------------------------------------
!
!  Close and compress velocity file
!
!-----------------------------------------------------------------------
!
  istatus=nf_close(ncid)
  CALL cmprsgz(cdfname)
END SUBROUTINE wtvrtiltcdf

SUBROUTINE wtvvtiltcdf(ngate,maxazim,nazim,fname,outdir,varname,       &
                       radname,radlat,radlon,radelv,ivcp,elv,          &
                       rmisval,rngfval,itimcdf,frtime,initime,         &
                       vnyquist,rfrgate,                               &
                       azim,beamw,gtspc,vvor)
!
! PURPOSE: Write one tilt of radar reflectivity data in NetCDF format.
! AUTHOR: Keith Brewster, CAPS
!
  IMPLICIT NONE
  INTEGER :: ngate
  INTEGER :: maxazim
  INTEGER :: nazim
  CHARACTER (LEN=256) :: fname
  CHARACTER (LEN=80)  :: outdir
  CHARACTER (LEN=40)  :: varname
  CHARACTER (LEN=4)   :: radname
  REAL :: radlat
  REAL :: radlon
  REAL :: radelv
  INTEGER :: ivcp
  REAL :: elv
  REAL :: rmisval
  REAL :: rngfval
  INTEGER :: itimcdf
  REAL :: frtime
  INTEGER :: initime
  REAL :: vnyquist
  REAL :: rfrgate
!
  REAL :: azim(maxazim)
  REAL :: beamw(maxazim)
  REAL :: gtspc(maxazim)
  REAL :: vvor(ngate,maxazim)
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
!  netCDF output variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: swpdim(2)
  INTEGER :: azmdim(1)
!
!-----------------------------------------------------------------------
!
!  Variable indexes and descriptors
!
!-----------------------------------------------------------------------
!
  INTEGER :: nazimid,ngateid
  INTEGER :: azmid,bmwid,elvid,gtwid,vorid
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: lrname,lvarnam,ifsecs
  INTEGER :: istr1d(1),iknt1d(1)
  INTEGER :: istr2d(2),iknt2d(2)
!
  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ncmode=nf_clobber
  lrname=LEN_TRIM(runname)
  lvarnam=LEN_TRIM(varname)
  istr1d(1)=1
  istr2d(1)=1
  istr2d(2)=1
  iknt1d(1)=nazim
  iknt2d(1)=ngate
  iknt2d(2)=nazim
!
!-----------------------------------------------------------------------
!
!  Open netCDF file
!
!-----------------------------------------------------------------------
!
  write(cdfname,'(a,a,a)') TRIM(outdir),'/',TRIM(fname)
  write(6,'(a,a)') ' Writing data to: ',TRIM(cdfname)
  istatus=nf_create(cdfname,ncmode,ncid)
!
!-----------------------------------------------------------------------
!
!  Define Dimensions
!
!-----------------------------------------------------------------------
!
  istatus=nf_def_dim(ncid,'Azimuth',nazim,nazimid)
  istatus=nf_def_dim(ncid,'Gate',ngate,ngateid)
  azmdim(1)=nazimid
  swpdim(1)=ngateid
  swpdim(2)=nazimid
!
!-----------------------------------------------------------------------
!
!  Define Variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_def_var(ncid,'Azimuth',nf_real,1,azmdim,azmid)
  istatus=nf_def_var(ncid,'BeamWidth',nf_real,1,azmdim,bmwid)
  istatus=nf_def_var(ncid,'GateWidth',nf_real,1,azmdim,gtwid)
  istatus=nf_def_var(ncid,TRIM(varname),nf_real,2,swpdim,vorid)
!
!-----------------------------------------------------------------------
!
!  Set Variable Attributes
!
!-----------------------------------------------------------------------
!
  istatus=nf_put_att_text(ncid,azmid,'Units',7,'Degrees')
  istatus=nf_put_att_text(ncid,bmwid,'Units',7,'Degrees')
  istatus=nf_put_att_text(ncid,gtwid,'Units',6,'Meters')
  istatus=nf_put_att_text(ncid,vorid,'Units',9,'1/Seconds')
!
!-----------------------------------------------------------------------
!
!  Set Global Attributes
!
!-----------------------------------------------------------------------
!
  istatus=nf_put_att_text(ncid,nf_global,'TypeName',lvarnam,TRIM(varname))
  istatus=nf_put_att_text(ncid,nf_global,'DataType',9,'RadialSet')
  istatus=nf_put_att_real(ncid,nf_global,'Latitude',nf_real,1,radlat)
  istatus=nf_put_att_real(ncid,nf_global,'Longitude',nf_real,1,radlon)
  istatus=nf_put_att_real(ncid,nf_global,'Height',nf_real,1,radelv)
  istatus=nf_put_att_int(ncid,nf_global,'Time',nf_int,1,itimcdf)
  istatus=nf_put_att_real(ncid,nf_global,'FractionalTime',nf_real,1,frtime)
  istatus=nf_put_att_text(ncid,nf_global,'attributes',26,  &
                            ' Nyquist_Vel radarName vcp')
  istatus=nf_put_att_text(ncid,nf_global,'Nyquist_Vel-unit',15,'MetersPerSecond')
  istatus=nf_put_att_real(ncid,nf_global,'Nyquist_Vel-value',nf_real,1,vnyquist)
  istatus=nf_put_att_text(ncid,nf_global,'radarName-unit',13,'dimensionless')
  istatus=nf_put_att_text(ncid,nf_global,'radarName-value',4,radname)
  istatus=nf_put_att_text(ncid,nf_global,'vcp-unit',13,'dimensionless')
  istatus=nf_put_att_int(ncid,nf_global,'vcp-value',nf_int,1,ivcp)
  istatus=nf_put_att_real(ncid,nf_global,'Elevation',nf_real,1,elv)
  istatus=nf_put_att_text(ncid,nf_global,'ElevationUnits',7,'Degrees')
  istatus=nf_put_att_real(ncid,nf_global,'RangeToFirstGate',nf_real,1,rfrgate)
  istatus=nf_put_att_text(ncid,nf_global,'RangeToFirstGateUnits',6,'Meters')
  istatus=nf_put_att_real(ncid,nf_global,'MissingData',nf_real,1,rmisval)
  istatus=nf_put_att_real(ncid,nf_global,'RangeFolded',nf_real,1,rngfval)

  istatus=nf_put_att_text(ncid,nf_global,'Runname',lrname,runname)
  istatus=nf_put_att_int(ncid,nf_global,'InitialTime',nf_int,1,initime)
  istatus=nf_put_att_real(ncid,nf_global,'ForecastSeconds',nf_real,1,curtim)
  
  istatus=nf_enddef(ncid)
!
!-----------------------------------------------------------------------
!
!  Write data variables
!
!-----------------------------------------------------------------------
!
  istatus=nf_put_vara_real(ncid,azmid,istr1d,iknt1d,azim)
  istatus=nf_put_vara_real(ncid,bmwid,istr1d,iknt1d,beamw)
  istatus=nf_put_vara_real(ncid,gtwid,istr1d,iknt1d,gtspc)
  istatus=nf_put_vara_real(ncid,vorid,istr2d,iknt2d,vvor)
!
!-----------------------------------------------------------------------
!
!  Close and compress reflectivity file
!
!-----------------------------------------------------------------------
!
  istatus=nf_close(ncid)
  CALL cmprsgz(cdfname)
END SUBROUTINE wtvvtiltcdf
