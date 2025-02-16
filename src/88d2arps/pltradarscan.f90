!
!##################################################################
!##################################################################
!######                                                      ######
!######      Advanced Regional Prediction System (ARPS)      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
PROGRAM PLT_RADAR_TILT
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!     Plot radar reflectivity and radial velocity in a tilt and
!              vertical cross-sections 
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  03/19/02.
!
!  MODIFICATION HISTORY:
!  06/12/2005 Keith Brewster, CAPS
!  Added reading of radar emulator volume files (radfmt=2)
!
!  11/15/2005 Keith Brewster, CAPS
!  Updated for changes to content of emulator volume files.
!
!  11/23/2005 Keith Brewster, CAPS
!  Restructured code to use considerably less memory.
!
!  10/19/2007 Keith Brewster, CAPS
!  Added max,min reflectivity to plotting NAMELIST.
!  Removed the grid namelist entirely as those variables as unused.
!  Fixed a bug im final deallocation logic due to uninitialized variable.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, PARAMETER :: mxradfile= 200
  INTEGER, PARAMETER :: mxmapfile= 10
  INTEGER, PARAMETER :: nref_tilts_max=16
  INTEGER, PARAMETER :: nvel_tilts_max=16
  INTEGER, PARAMETER :: n_elevation_max=16
  INTEGER, PARAMETER :: n_azimuth_max=800

!
!-----------------------------------------------------------------------
!
!  Namelist 
!
!-----------------------------------------------------------------------
!

  CHARACTER (LEN=256) :: filepath
  INTEGER :: radfmt
  INTEGER :: numfiles
  CHARACTER (LEN=256) :: files(mxradfile)
  NAMELIST /input_data/ filepath,radfmt,numfiles,files

  CHARACTER(LEN=4) :: radar_symbol
  NAMELIST /radar_index/ radar_symbol

  INTEGER :: mapproj
  REAL :: trulat1
  REAL :: trulat2
  REAL :: trulon
  REAL :: sclfct
  NAMELIST /projection/ mapproj,trulat1,trulat2,trulon,sclfct

  INTEGER :: overlayrefcontour
  INTEGER :: pltref,pltzdr,pltkdp,pltrhohv
  INTEGER :: pltqr,pltqs,pltqh
  INTEGER :: pltu,pltv,pltw,pltt
  INTEGER :: n_elevation_ref
  INTEGER :: n_elevation_vel
  INTEGER :: n_azimuth_ref
  INTEGER :: n_azimuth_vel
  REAL :: elevations_ref(n_elevation_max)
  REAL :: elevations_vel(n_elevation_max)
  REAL :: azimuth_ref(n_azimuth_max)
  REAL :: azimuth_vel(n_azimuth_max)
  REAL :: h_range_max,h_range_max_xs
  REAL :: h_range_max_e,h_range_max_w,h_range_max_s,h_range_max_n
  REAL :: z_range_max
  REAL :: elev_error_range
  REAL :: pltmin_ref
  REAL :: pltmax_ref
  INTEGER :: coltab_ref
  INTEGER :: ibgncol_ref
  INTEGER :: iendcol_ref
  INTEGER :: rhimode
  CHARACTER(LEN=256) :: coltabfn_ref
  REAL :: pltmin_vel
  REAL :: pltmax_vel
  INTEGER :: coltab_vel
  INTEGER :: ibgncol_vel
  INTEGER :: iendcol_vel
  CHARACTER(LEN=256) :: coltabfn_vel
  INTEGER :: coltab_dpv
  INTEGER :: ibgncol_dpv
  INTEGER :: iendcol_dpv
  CHARACTER(LEN=256) :: coltabfn_dpv
  NAMELIST /plotting/ h_range_max_e,h_range_max_w,h_range_max_s,        &
                      h_range_max_n,z_range_max,h_range_max_xs,         &
                      elev_error_range,overlayrefcontour,rhimode,       &
                      pltref,pltzdr,pltkdp,pltrhohv,                    &
                      pltqr,pltqs,pltqh,                                &
                      pltu,pltv,pltw,pltt,                              &
                      n_elevation_ref,elevations_ref,                   &
                      n_elevation_vel,elevations_vel,                   &
                      n_azimuth_ref,azimuth_ref,                        &
                      n_azimuth_vel,azimuth_vel,                        &
                      pltmin_ref,pltmax_ref,                            &
                      coltab_ref,ibgncol_ref,iendcol_ref,coltabfn_ref,  &
                      pltmin_vel,pltmax_vel,                            &
                      coltab_vel,ibgncol_vel,iendcol_vel,coltabfn_vel,  &
                      coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv

  INTEGER :: ovrmap
  INTEGER :: mapgrid
  INTEGER :: nmapfile
  INTEGER :: mapcol(mxmapfile)
  INTEGER :: mapline_style(mxmapfile)
  INTEGER :: mapgridcol
  REAL :: latgrid
  REAL :: longrid
  CHARACTER (LEN=256) :: mapfile(mxmapfile)
  NAMELIST /map_plot/ ovrmap,mapgrid,latgrid,longrid,mapgridcol, &
                      nmapfile,mapfile,mapcol,mapline_style
  INTEGER :: lmapfile
  REAL :: xorig, yorig
!
!-----------------------------------------------------------------------
!
! plot variables
!
!-----------------------------------------------------------------------
!
  REAL :: cl(100)

  REAL, ALLOCATABLE :: xrplt(:,:)
  REAL, ALLOCATABLE :: yrplt(:,:)
  REAL, ALLOCATABLE :: xvplt(:,:)
  REAL, ALLOCATABLE :: yvplt(:,:)
  REAL, ALLOCATABLE :: rplt(:,:)
  REAL, ALLOCATABLE :: rpltm(:,:)
  REAL, ALLOCATABLE :: zplt(:,:)
  REAL, ALLOCATABLE :: xrw(:)
  REAL, ALLOCATABLE :: yrw(:)
  REAL, ALLOCATABLE :: xvw(:)
  REAL, ALLOCATABLE :: yvw(:)
  REAL, ALLOCATABLE :: elevmean(:)

  INTEGER, ALLOCATABLE :: irwxy(:,:)
  INTEGER, ALLOCATABLE :: ivwxy(:,:)
  INTEGER, ALLOCATABLE :: iwrz(:,:)

  INTEGER :: ncl
  INTEGER :: mode
!
!-----------------------------------------------------------------------
!
!  INPUT RADAR PARAMETER AND OBSERVATION 
!
!-----------------------------------------------------------------------
!
  REAL :: radar_lat
  REAL :: radar_lon
  REAL :: radar_elevation

  INTEGER :: iyear,imon,iday,ihour,imin,isec

  CHARACTER(LEN=256) :: datafile
  CHARACTER(LEN=80)  :: runname
  INTEGER :: lfnkey,abstsec
  INTEGER :: wrtuaref,wrtuavel
  INTEGER :: wrtvort,wrtdualp,wrtqx,wrtuvwt
  INTEGER :: idummy,vcp
  REAL :: beamwid,rmisval,rngfval,curtim

  REAL, ALLOCATABLE :: elvref(:)
  REAL, ALLOCATABLE :: elvvel(:)
  REAL, ALLOCATABLE :: reftilt(:,:)
  REAL, ALLOCATABLE :: vartilt(:,:)
  REAL, ALLOCATABLE :: veltilt(:,:)
  REAL, ALLOCATABLE :: datarz(:,:)

  INTEGER :: maxvelgate
  INTEGER :: maxrefgate
  INTEGER :: maxazim
  INTEGER :: maxelev

  INTEGER, ALLOCATABLE :: itimvol(:)

  REAL, ALLOCATABLE :: rngvvol(:,:)
  REAL, ALLOCATABLE :: azmvvol(:,:)
  REAL, ALLOCATABLE :: elvvvol(:,:)
  REAL, ALLOCATABLE :: velvol(:,:,:)

  REAL, ALLOCATABLE :: rngrvol(:,:)
  REAL, ALLOCATABLE :: azmrvol(:,:)
  REAL, ALLOCATABLE :: elvrvol(:,:)
  REAL, ALLOCATABLE :: refvol(:,:,:)
  REAL, ALLOCATABLE :: vortvol(:,:,:)

  REAL, ALLOCATABLE :: zdrvol(:,:,:)
  REAL, ALLOCATABLE :: kdpvol(:,:,:)
  REAL, ALLOCATABLE :: rhohvvol(:,:,:)
  REAL, ALLOCATABLE :: hailflg(:,:,:)
  REAL, ALLOCATABLE :: qrvol(:,:,:)
  REAL, ALLOCATABLE :: qsvol(:,:,:)
  REAL, ALLOCATABLE :: qhvol(:,:,:)
  REAL, ALLOCATABLE :: uvol(:,:,:)
  REAL, ALLOCATABLE :: vvol(:,:,:)
  REAL, ALLOCATABLE :: wvol(:,:,:)
  REAL, ALLOCATABLE :: tvol(:,:,:)

  REAL, ALLOCATABLE :: vnyquist(:)
  REAL, ALLOCATABLE :: uavol(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: ifiles   
  CHARACTER(LEN=30) :: varname
  CHARACTER(LEN=100) :: string 
  CHARACTER(LEN=100) :: stringtime
  REAL :: radmax ,rad
  INTEGER :: ialstatus
  INTEGER :: i_angle,i_tilt,igate
  INTEGER :: j_razim_start,nref_tilts,nvel_tilts
  INTEGER :: iazimuth
  REAL :: razim_start
  INTEGER :: gsp_vel,rfrst_vel

  REAL :: elev_to_plot, xscale, yscale
  REAL :: height,sfcrng

  INTEGER :: kref,kvel,ktilt
  INTEGER :: kref_tilts,kvel_tilts
  REAL :: refmax, refmin, reflmin, reflmax
  REAL :: velmax, velmin
  REAL :: eleva, deg2rad

  REAL :: xradar, yradar, tem0,tem1,tem2,tem3,tem4,xlab
  LOGICAL :: fexist
  INTEGER :: zxplot_initialized = 0
  REAL :: missing_value = -999.0
  INTEGER :: iunit

  INTEGER :: i,j,k,jj

  CHARACTER(LEN=256) :: psfilename
  INTEGER :: ilen,ialloc

  REAL :: factorx,factory
  REAL :: delazm,dazmin,xtick,ztick,zradar
  REAL :: xmin,xmax,zmin,zmax
  REAL :: elvsum,data_azim,last_azim,azmin,xconst,yconst
  INTEGER :: elvknt,jidx,jazim,nazim,rhidata

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  deg2rad=acos(-1.0)/180.0

  sclfct = 1.0
  rhimode = 0
  wrtuaref = 0
  wrtuavel = 0
  wrtvort  = 0
  wrtdualp = 0
  wrtqx = 0
  wrtuvwt  = 0

  pltmin_ref=5.
  pltmax_ref=80.
  pltmin_vel=-150.
  pltmax_vel=150.
  coltab_ref  = 3
  ibgncol_ref = 32
  iendcol_ref = 46
  coltabfn_ref='./data/arpsplt/hubcaps5.pltcbar'
  coltab_vel  = -1
  ibgncol_vel = 2
  iendcol_vel = 21
  coltabfn_vel='./data/arpsplt/velocity_solo_real.table'
  coltab_dpv  = 1
  ibgncol_dpv = 3
  iendcol_dpv = 24
  coltabfn_dpv='./data/arpsplt/hubcaps5.pltcbar'
!
!-----------------------------------------------------------------------
! Input data file name and grid 
!-----------------------------------------------------------------------

  READ(5,input_data,ERR=100)
  WRITE(6,'(a)') 'Namelist input_data was sucessfully read in.'

  READ(5,radar_index,ERR=100) 
  WRITE(6,'(a)') 'Namelist radar_index was sucessfully read in.'

!-----------------------------------------------------------------------
! Input map projection parameters
!-----------------------------------------------------------------------
 
  READ (5,projection,ERR=100)
  WRITE(6,'(a)') 'Namelist projection was sucessfully read in.'

  READ(5,plotting,ERR=100)
  WRITE(6,'(a)') 'Namelist plotting sucessfully read in.'

  READ(5,map_plot,ERR=100)
  WRITE(6,'(a)') 'Namelist map_plot sucessfully read in.'

!-----------------------------------------------------------------------
! Read in radar data volume 
!-----------------------------------------------------------------------

  DO ifiles=1,numfiles
!
!  read in data on a tilt
!
    WRITE(datafile,'(a,a,a)') TRIM(filepath),'/',TRIM(files(ifiles))
    INQUIRE(FILE=TRIM(datafile),EXIST = fexist)
    IF( .not. fexist) THEN
      WRITE(6,'(a,a,a)') 'Radar data file ',TRIM(datafile),' not found.'
      CYCLE
    END IF

    WRITE(6,'(a,i4)') ' Radar format = ',radfmt
    iunit=27
    IF( radfmt == 1 ) THEN
      WRITE(6,'(a,a)') ' Opening: ',TRIM(datafile)
      OPEN(iunit,file=TRIM(datafile),form='unformatted',status='old')

      READ(iunit) maxvelgate,maxrefgate,maxazim,maxelev
      READ(iunit) radar_elevation,radar_lat, radar_lon
      READ(iunit) iyear,imon,iday,ihour,imin,isec

      WRITE(6,'(3i10)') maxvelgate,maxrefgate,maxazim,maxelev 
      WRITE(6,'(3f10.2)') radar_elevation,radar_lat,radar_lon
      WRITE(6,'(6i6)') iyear,imon,iday,ihour,imin,isec

      WRITE(6,'(a)') 'Allocating radar data space'

      ALLOCATE( rngrvol(maxrefgate,maxelev), STAT=ialloc )
      IF(ialloc /= 0) THEN
        WRITE(6,'(a,i6)') ' Error allocating rngrvol: ',ialloc
        STOP
      END IF
      ALLOCATE( azmrvol(maxazim,maxelev), STAT=ialloc )
      IF(ialloc /= 0) THEN
        WRITE(6,'(a,i6)') ' Error allocating azmrvol: ',ialloc
        STOP
      END IF
      ALLOCATE( elvrvol(maxazim,maxelev), STAT=ialloc )
      IF(ialloc /= 0) THEN
        WRITE(6,'(a,i6)') ' Error allocating elvrvol: ',ialloc
        STOP
      END IF
      ALLOCATE( refvol(maxrefgate,maxazim,maxelev), STAT=ialloc )
      IF(ialloc /= 0) THEN
        WRITE(6,'(a,i6)') ' Error allocating refvol: ',ialloc
        STOP
      END IF

      ALLOCATE( rngvvol(maxvelgate,maxelev), STAT=ialloc )
      IF(ialloc /= 0) THEN
        WRITE(6,'(a,i6)') ' Error allocating rngvvol: ',ialloc
        STOP
      END IF
      ALLOCATE( azmvvol(maxazim,maxelev), STAT=ialloc )
      IF(ialloc /= 0) THEN
        WRITE(6,'(a,i6)') ' Error allocating azmvvol: ',ialloc
        STOP
      END IF
      ALLOCATE( elvvvol(maxazim,maxelev), STAT=ialloc )
      IF(ialloc /= 0) THEN
        WRITE(6,'(a,i6)') ' Error allocating elvvvol: ',ialloc
        STOP
      END IF
      ALLOCATE( velvol(maxvelgate,maxazim,maxelev), STAT=ialloc )
      IF(ialloc /= 0) THEN
        WRITE(6,'(a,i6)') ' Error allocating velvol: ',ialloc
        STOP
      END IF
      WRITE(6,'(a)') 'Success allocating radar data space'

      READ(iunit) rngrvol
      READ(iunit) azmrvol
      READ(iunit) elvrvol
      READ(iunit) refvol
  
      READ(iunit) rngvvol
      READ(iunit) azmvvol
      READ(iunit) elvvvol
      READ(iunit) velvol

      close(iunit)

      psfilename=TRIM(files(ifiles))
      ilen=LEN(TRIM(psfilename))
      iday=ICHAR(psfilename(ilen-12:ilen-12))-ICHAR("0")
      iday=iday*10+(ICHAR(psfilename(ilen-11:ilen-11))-ICHAR("0"))
      ihour=ICHAR(psfilename(ilen-9:ilen-9))-ICHAR("0")
      ihour=ihour*10+(ICHAR(psfilename(ilen-8:ilen-8))-ICHAR("0"))
      imin=ICHAR(psfilename(ilen-7:ilen-7))-ICHAR("0")
      imin=imin*10+(ICHAR(psfilename(ilen-6:ilen-6))-ICHAR("0"))

    ELSE IF (radfmt == 2) THEN   ! radar emulator binary file
      WRITE(6,'(a,a)') ' Opening: ',TRIM(datafile)
      OPEN(iunit,FILE=TRIM(datafile),FORM='unformatted',STATUS='old')
      READ(iunit) runname
      lfnkey=80
      CALL gtlfnkey(runname,lfnkey)
      WRITE(6,'(a,a)') ' Runname:',runname(1:lfnkey)
      READ(iunit) radar_symbol
      WRITE(6,'(a,a)') ' Radar name:',radar_symbol
      READ(iunit) maxvelgate,maxrefgate,maxazim,maxelev
      WRITE(6,'(a,i6,a,i6,a,i6)')                                          &
        ' Ngate : ',maxvelgate,' Nazim: ',maxazim,' maxelev: ',maxelev
      READ(iunit) iyear,imon,iday,ihour,imin,isec
      WRITE(6,'(2(a,i2.2),a,i4.4,3(a,i2.2))')                              &
        ' Date: ',imon,'/',iday,'/',iyear,' Time: ',ihour,':',imin,':',isec
      READ(iunit) radar_elevation,radar_lat,radar_lon
      WRITE(6,'(a,f9.2,a,f9.2,a,f9.2)')                                      &
        ' Latitude:',radar_lat,'  Longitude:',radar_lon, &
        ' Elevation:',radar_elevation
      READ(iunit) vcp,beamwid,rmisval,rngfval,curtim
      WRITE(6,'(a,i5,a,f9.2,a,f10.1,a,f10.1)')                               &
        ' VCP:',vcp,' Beamwidth:',beamwid,                                   &
        ' MissVal: ',rmisval,' RngfVal :',rngfval
!
      CALL ctim2abss(iyear,imon,iday,ihour,imin,isec,abstsec)
      READ(iunit) wrtuaref,wrtuavel,wrtvort,wrtdualp,                        &
                  wrtqx,wrtuvwt,idummy,idummy
      WRITE(6,'(a,i4)') ' wrtuaref =',wrtuaref
      WRITE(6,'(a,i4)') ' wrtuavel =',wrtuavel
      WRITE(6,'(a,i4)') ' wrtvort  =',wrtvort
      WRITE(6,'(a,i4)') ' wrtdualp =',wrtdualp
      WRITE(6,'(a,i4)') ' wrtqx    =',wrtqx
      WRITE(6,'(a,i4)') ' wrtuvwt  =',wrtuvwt

      WRITE(6,'(a)') ' Allocating volume memory '

      ALLOCATE(itimvol(maxelev),STAT=ialloc)
      ALLOCATE(rngrvol(maxrefgate,maxelev),STAT=ialloc)
      ALLOCATE(azmrvol(maxazim,maxelev),STAT=ialloc)
      ALLOCATE(elvrvol(maxazim,maxelev),STAT=ialloc)
      ALLOCATE(refvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
      IF(ialloc /= 0) THEN
        WRITE(6,'(a,i6)') 'Error allocating refvol:',ialloc
        STOP
      END IF
      IF(wrtuaref /= 0 .OR. wrtuavel /= 0 ) &
        ALLOCATE(uavol(maxrefgate,maxazim,maxelev), STAT=ialloc)
      IF(ialloc /= 0) THEN
        WRITE(6,'(a,i6)') 'Error allocating uavol:',ialloc
        STOP
      END IF

      IF(wrtdualp /= 0) THEN
        ALLOCATE(zdrvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
        ALLOCATE(kdpvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
        ALLOCATE(rhohvvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
        ALLOCATE(hailflg(maxrefgate,maxazim,maxelev),STAT=ialloc)
        IF(ialloc /= 0) THEN
          WRITE(6,'(a,i6)') 'Error allocating dual-pol vol:',ialloc
          STOP
        END IF
      END IF

      IF(wrtqx /= 0) THEN
        ALLOCATE(qrvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
        ALLOCATE(qsvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
        ALLOCATE(qhvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
        IF(ialloc /= 0) THEN
          WRITE(6,'(a,i6)') 'Error allocating qx vol:',ialloc
          STOP
        END IF
      END IF

      IF(wrtuvwt /= 0) THEN
        ALLOCATE(uvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
        ALLOCATE(vvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
        ALLOCATE(wvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
        ALLOCATE(tvol(maxrefgate,maxazim,maxelev),STAT=ialloc)
        IF(ialloc /= 0) THEN
          WRITE(6,'(a,i6)') 'Error allocating uvwt vol:',ialloc
          STOP
        END IF
      END IF

      ALLOCATE(vnyquist(maxelev))
      ALLOCATE(rngvvol(maxvelgate,maxelev))
      ALLOCATE(azmvvol(maxazim,maxelev))
      ALLOCATE(elvvvol(maxazim,maxelev))
      ALLOCATE(velvol(maxvelgate,maxazim,maxelev),STAT=ialloc)
      IF(ialloc == 0) WRITE(6,'(a)') ' Vel volume memory allocated OK'

      IF(wrtvort /= 0) THEN
        ALLOCATE(vortvol(maxvelgate,maxazim,maxelev),STAT=ialloc)
        IF(ialloc == 0) WRITE(6,'(a)') ' Vort volume memory allocated OK'
      END IF
                     
      WRITE(6,'(a)') ' Reading Time Array'
      READ(iunit) itimvol
      WRITE(6,'(a)') ' Reading Reflectivity Arrays'
      READ(iunit) rngrvol
      READ(iunit) azmrvol
      READ(iunit) elvrvol
      READ(iunit) refvol
      IF( wrtuaref /= 0 ) READ(iunit) uavol
      IF( wrtdualp /= 0 ) THEN
        WRITE(6,'(a)') ' Reading Dual-Pol Arrays'
        READ(iunit) zdrvol
        READ(iunit) kdpvol
        READ(iunit) rhohvvol
      END IF

      WRITE(6,'(a)') ' Reading Velocity Arrays'
      READ(iunit) vnyquist
      READ(iunit) rngvvol
      READ(iunit) azmvvol
      READ(iunit) elvvvol
      READ(iunit) velvol
      IF( wrtuavel /= 0 ) READ(iunit) uavol
      IF( wrtvort /= 0 ) READ(iunit) vortvol

      IF( wrtqx /= 0 ) THEN
        WRITE(6,'(a)') ' Reading qr/qs/qh Arrays'
        READ(iunit) qrvol
        READ(iunit) qsvol
        READ(iunit) qhvol
      END IF

      IF( wrtuvwt /= 0 ) THEN
        WRITE(6,'(a)') ' Reading u,v,w,t Arrays'
        READ(iunit) uvol
        READ(iunit) vvol
        READ(iunit) wvol
        READ(iunit) tvol
      END IF

!
      WRITE(6,'(a)') ' Reading successfully completed'
      CLOSE(iunit)

    ELSE IF (radfmt == 3) THEN   ! Radar unfolding diagnostic file

      maxrefgate = 1
      maxelev = 1
      OPEN(iunit,file=TRIM(datafile),form='unformatted',status='old')
      READ(iunit) radar_symbol
      READ(iunit) iyear,imon,iday,ihour,imin,isec

      READ(iunit) i_angle,i_tilt,rfrst_vel,gsp_vel
      READ(iunit) maxvelgate,maxazim
      READ(iunit) radar_lat,radar_lon,radar_elevation
      WRITE(6,'(a,f9.2,a,f9.2,a,f9.2)')                  &
        ' Latitude:',radar_lat,'  Longitude:',radar_lon, &
        ' Elevation:',radar_elevation

      WRITE(6,'(a)') 'Allocating radar data space'

      ALLOCATE(rngrvol(maxrefgate,1))
      ALLOCATE(azmrvol(maxazim,1))
      ALLOCATE(elvrvol(maxazim,1))
      ALLOCATE(refvol(maxrefgate,maxazim,1))
      refvol=-999.

      ALLOCATE(vnyquist(maxelev))
      ALLOCATE(rngvvol(maxvelgate,maxelev))
      ALLOCATE(azmvvol(maxazim,maxelev))
      ALLOCATE(elvvvol(maxazim,maxelev))
      ALLOCATE(velvol(maxvelgate,maxazim,maxelev),STAT=ialloc)
      IF(ialloc == 0) WRITE(6,'(a)') ' Vel volume memory allocated OK'

      READ(iunit) vnyquist
      READ(iunit) elvvvol
      READ(iunit) azmvvol
      elvrvol=elvvvol
      azmrvol=azmvvol

      READ(iunit) velvol

      DO igate=1,maxvelgate
        rngvvol(igate,1)=rfrst_vel+(igate-1)*gsp_vel
      END DO
      rngrvol=rfrst_vel

    ELSE IF (radfmt == 4) THEN   ! CASA-formatted NetCDF file
      WRITE(6,'(a,i6,a)') ' Radar format number ',radfmt,' not yet implemented'
      STOP
      
    ELSE
      WRITE(6,'(a,i6,a)') ' Radar format number ',radfmt,' not supported.'
      STOP
    END IF

    IF ( mapproj == 0 ) THEN
        trulat1 = radar_lat
        trulat2 = radar_lat
        trulon  = radar_lon
    END IF

!-----------------------------------------------------------------------
! Data variable unit conversions
!-----------------------------------------------------------------------

   IF( wrtqx > 0) THEN
     DO k=1,maxelev
       DO j=1,maxazim
         DO i=1,maxrefgate
           qrvol(i,j,k)=1000.*qrvol(i,j,k)
           qsvol(i,j,k)=1000.*qsvol(i,j,k)
           qhvol(i,j,k)=1000.*qhvol(i,j,k)
         END DO
       END DO
     END DO
   END IF

   IF( wrtuvwt > 0) THEN
     DO k=1,maxelev
       DO j=1,maxazim
         DO i=1,maxrefgate
           tvol(i,j,k)=tvol(i,j,k)-273.15
         END DO
       END DO
     END DO
   END IF

   IF(wrtdualp /= 0) THEN
     hailflg=0.
     DO k=1,maxelev
       DO j=1,maxazim
         DO i=1,maxrefgate
           IF (refvol(i,j,k) >= 45. ) THEN
             IF( zdrvol(i,j,k) < 2.0) THEN
               hailflg(i,j,k)=1.0
             ELSE IF(zdrvol(i,j,k) < 2.25) THEN
               hailflg(i,j,k)=0.95
             ELSE IF(zdrvol(i,j,k) < 2.5) THEN
               hailflg(i,j,k)=0.90
             ELSE IF(zdrvol(i,j,k) < 2.75) THEN
               hailflg(i,j,k)=0.80
             ELSE IF(zdrvol(i,j,k) < 3.0) THEN
               hailflg(i,j,k)=0.60
             ELSE IF(zdrvol(i,j,k) < 3.25) THEN
               hailflg(i,j,k)=0.50
             ELSE
               hailflg(i,j,k)=0.2
             END IF
           END IF
         END DO
       END DO
     END DO
   END IF
     
!-----------------------------------------------------------------------
! plot set-up
!-----------------------------------------------------------------------

    rhidata=0
    IF(vcp == -99) THEN
      WRITE(6,'(a,i,/a)') '  Data are RHI mode, vcp =',vcp,&
                          '  Skipping PPI plots'
      rhidata=1
      rhimode=1
    END IF
    ALLOCATE(elvref(maxelev))
    ALLOCATE(elvvel(maxelev))

    factorx=0.3
    factory=0.3
    IF(h_range_max_w + h_range_max_e > h_range_max_s + h_range_max_n ) THEN
       factory=factorx*(h_range_max_s + h_range_max_n)/(h_range_max_w + h_range_max_e)
       h_range_max=(h_range_max_e+h_range_max_w)/2
    ELSE
       factorx=factory*(h_range_max_w + h_range_max_e)/(h_range_max_s + h_range_max_n)
       h_range_max=(h_range_max_n+h_range_max_s)/2
    END IF

    IF( n_elevation_ref > 0 .or. n_azimuth_ref > 0 .or.                    &
        n_elevation_vel > 0 .or. n_azimuth_vel > 0 ) THEN
      IF( zxplot_initialized == 0 ) THEN 
        CALL xpsfn('zxout'//'.ps', 5) 
                    ! PS output will be called zxout.ps
        CALL xpaprlnth( 1.0 )

        CALL xdevic     ! begin graphic
        zxplot_initialized = 1 
      END IF
    END IF

    elvref=-99.
    DO ktilt=1,maxelev
      elvsum=0.
      elvknt=0
      DO j=1,maxazim
        IF(elvrvol(j,ktilt) > -91.) THEN
          elvknt=elvknt+1
          elvsum=elvsum+elvrvol(j,ktilt)
        END IF
      END DO
      IF(elvknt > 0) THEN
        elvref(ktilt)=elvsum/float(elvknt)
        nref_tilts=ktilt
      ELSE
        elvref(ktilt)=-99.
      END IF
      WRITE(6,'(a,i4,a,f10.2)') ' Elevation(',ktilt,') = ',elvref(ktilt)
    END DO
    WRITE(6,'(a,i4)') ' nref_tilts=',nref_tilts

    elvvel=-99.
    DO ktilt=1,maxelev
      elvsum=0.
      elvknt=0
      DO j=1,maxazim
        IF(elvvvol(j,ktilt) > -9.) THEN
          elvknt=elvknt+1
          elvsum=elvsum+elvvvol(j,ktilt)
        END IF
      END DO
      IF(elvknt > 0) THEN
        elvvel(ktilt)=elvsum/float(elvknt)
        nvel_tilts=ktilt
      ELSE
        elvvel(ktilt)=-99.
      END IF
      WRITE(6,'(a,i4,a,f10.2)') ' Elevation(',ktilt,') = ',elvvel(ktilt)
    END DO
    WRITE(6,'(a,i4)') ' nvel_tilts=',nvel_tilts

!-----------------------------------------------------------------------
! plot reflectivity variables
!-----------------------------------------------------------------------

    IF(rhidata < 1 .AND. n_elevation_ref > 0) THEN

    WRITE(6,'(a)') ' Allocating xrplt,yrplt,reftilt ...'
    ALLOCATE( xrplt(maxrefgate,(maxazim+1)) )
    ALLOCATE( yrplt(maxrefgate,(maxazim+1)) )
    ALLOCATE( reftilt(maxrefgate,(maxazim+1)) )
    ALLOCATE( vartilt(maxrefgate,(maxazim+1)) )
    ALLOCATE( irwxy(maxrefgate,(maxazim+1)) )
    ALLOCATE( xrw(8*maxrefgate) )
    ALLOCATE( yrw(8*maxrefgate) )

    WRITE(6,'(a)') ' Reflectivity tilt plotting...'

    DO ktilt = 1,n_elevation_ref

      elev_to_plot = elevations_ref(ktilt)

      kref = 0 
      DO k = 1,maxelev
        IF( abs( elvref(k)-elev_to_plot )                   &
                            < elev_error_range ) THEN
          kref = k
          EXIT
        END IF
      END DO

      IF(kref == 0 ) CYCLE  ! no elevation angle 
                            ! within error range was found
      refmin = 1000.0
      refmax =-1000.0
      velmin = 1000.0
      refmax =-1000.0

      eleva=elvref(kref)

      IF(overlayrefcontour > 0 .OR. pltref > 0) THEN
        last_azim=-999.
        jazim=0
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            reflmax=-999.
            reflmin=999.
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              reftilt(i,jazim)=max(refvol(i,jidx,kref),-30.)

              IF(refvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(refvol(i,jidx,kref),refmax)
                refmin=min(refvol(i,jidx,kref),refmin)
                reflmax=max(refvol(i,jidx,kref),reflmax)
                reflmin=min(refvol(i,jidx,kref),reflmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        WRITE(6,'(a,2f10.2)') ' Refmin, Refmax on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        IF(pltref > 0) THEN
          CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
          CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

          CALL xlabmask(1)

          IF( coltab_ref == -1 ) CALL xstctfn(coltabfn_ref)
          CALL xsetclrs(coltab_ref)
          CALL xcolor(1)    ! Set color to black
          CALL xcfont(3)
!
          call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
          cl(1)=0.0
          cl(2)=5.0           ! Set tentative contour interval
          CALL xnctrs(2,20)  ! Set lower and upper limits
                              ! on the number of contours allowed
                              ! It's used by both xcolfil and xclimit
          mode=2

          call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
          call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
          CALL xctrclr(ibgncol_ref,iendcol_ref)   ! Use colors between 10 and 39 inclusive

          CALL xctrlim(pltmin_ref,pltmax_ref) ! Only contours above zero are plotted

          CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
          WRITE(6,'(a,f10.2)') &
       ' Plotting colfil reflectivity field at elevation angle = ', eleva
          CALL xcolfil(reftilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

          IF ( overlayrefcontour == 1 ) THEN
            call XCLTYP(0)
            CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
            mode=3
            ncl=1
            cl(1)=15.0
            CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                      maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
            cl(1)=30.0
            CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                      maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
            cl(1)=45.0
            CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                      maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
            cl(1)=55.0
            CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                      maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          END IF

          CALL xwdwof
          CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                      -h_range_max_s, h_range_max_n, 0.1*h_range_max)

          CALL xchsiz(0.025*2*h_range_max )  ! Set character size
          CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

          WRITE(string,'(a,f5.2)') 'Reflectivity at Elevation Angle=',eleva
          CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.04, trim(string) )

          WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
                imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
          CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

          CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
          CALL xcpalet(1)            ! Plot horizontal color palette

          xscale = h_range_max
          yscale = h_range_max
          xorig = -(h_range_max_e+h_range_max_w)/2.0
          yorig = -(h_range_max_s+h_range_max_n)/2.0

          xradar = 0.0
          yradar = 0.0
          CALL xthick(3)
          CALL xpenup( xradar-xscale*0.02, yradar )
          CALL xpendn( xradar+xscale*0.02, yradar )
          CALL xpenup( xradar, yradar-yscale*0.02 )
          CALL xpendn( xradar, yradar+yscale*0.02 )
          call xthick(1)
          CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
          CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

          xscale = h_range_max_e+h_range_max_w
          yscale = h_range_max_s+h_range_max_n
          CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                        xscale,yscale,xorig,yorig)
          CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

          IF( ovrmap /= 0 ) THEN
            DO i=1,MIN(mxmapfile,nmapfile)
              lmapfile=len_trim(mapfile(i))
              WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

              INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
              IF( .NOT.fexist) THEN
                WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                   ' not found. Corresponding map no plotted.'
              ELSE
                CALL xcolor(mapcol(i))
                IF(mapline_style(i) == 1) THEN
                  CALL xthick(1)
                  CALL xbrokn(6,3,6,3)
                ELSE IF(mapline_style(i) == 2) THEN
                  CALL xthick(1)
                ELSE IF(mapline_style(i) == 3) THEN
                  CALL xthick(3)
                  CALL xfull
                END IF
                CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
              END IF
            END DO
          END IF
    
          CALL xwdwof
          CALL xframe
        END IF  ! pltref
      END IF  ! pltref OR overlayref

      IF(wrtdualp > 0 .AND. pltzdr > 0) THEN
        last_azim=-999.
        jazim=0
        refmin = 1000.0
        refmax =-1000.0
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              vartilt(i,jazim)=min(max(zdrvol(i,jidx,kref),-10.),20.)

              IF(zdrvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(zdrvol(i,jidx,kref),refmax)
                refmin=min(zdrvol(i,jidx,kref),refmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        WRITE(6,'(a,2f10.2)') ' Zdrmin, Zdrmax on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
        CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        CALL xlabmask(1)

        IF( coltab_dpv == -1 ) CALL xstctfn(coltabfn_dpv)
        CALL xsetclrs(coltab_dpv)
        CALL xcolor(1)    ! Set color to black
        CALL xcfont(3)
!
        call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
        cl(1)=0.0
        cl(2)=0.5           ! Set tentative contour interval
        CALL xnctrs(5,50)   ! Set lower and upper limits
        mode=2

        call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
        call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
        CALL xctrclr(ibgncol_dpv,iendcol_dpv)

        CALL xctrlim(-2.0,6.0)

        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
        WRITE(6,'(a,f10.2)') &
       ' Plotting colfil Zdr field at elevation angle = ', eleva
        CALL xcolfil(vartilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

        IF ( overlayrefcontour == 1 ) THEN
          call XCLTYP(0)
          CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF

        CALL xwdwof
        CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

        CALL xchsiz(0.025*2*h_range_max )  ! Set character size
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

        WRITE(string,'(a,f5.2)') 'Zdr at Elevation Angle=',eleva
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
             h_range_max_n+2*h_range_max*0.04, trim(string) )

        WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcpalet(1)            ! Plot horizontal color palette

        xscale = h_range_max
        yscale = h_range_max
        xorig = -(h_range_max_e+h_range_max_w)/2.0
        yorig = -(h_range_max_s+h_range_max_n)/2.0

        xradar = 0.0
        yradar = 0.0
        CALL xthick(3)
        CALL xpenup( xradar-xscale*0.02, yradar )
        CALL xpendn( xradar+xscale*0.02, yradar )
        CALL xpenup( xradar, yradar-yscale*0.02 )
        CALL xpendn( xradar, yradar+yscale*0.02 )
        call xthick(1)
        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

        xscale = h_range_max_e+h_range_max_w
        yscale = h_range_max_s+h_range_max_n
        CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                    xscale,yscale,xorig,yorig)
        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        IF( ovrmap /= 0 ) THEN
          DO i=1,MIN(mxmapfile,nmapfile)
            lmapfile=len_trim(mapfile(i))
            WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

            INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
            IF( .NOT.fexist) THEN
              WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
            ELSE
              CALL xcolor(mapcol(i))
              IF(mapline_style(i) == 1) THEN
                CALL xthick(1)
                CALL xbrokn(6,3,6,3)
              ELSE IF(mapline_style(i) == 2) THEN
                CALL xthick(1)
              ELSE IF(mapline_style(i) == 3) THEN
                CALL xthick(3)
                CALL xfull
              END IF
              CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
            END IF
          END DO
        END IF
    
        CALL xwdwof
        CALL xframe
      END IF  ! pltzdr

      IF(wrtdualp > 0 .AND. pltkdp > 0) THEN
        last_azim=-999.
        jazim=0
        refmin = 1000.0
        refmax =-1000.0
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              vartilt(i,jazim)=max(kdpvol(i,jidx,kref),-30.)

              IF(kdpvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(kdpvol(i,jidx,kref),refmax)
                refmin=min(kdpvol(i,jidx,kref),refmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        WRITE(6,'(a,2f10.2)') ' Kdpmin, Kdpmax on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
        CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        CALL xlabmask(1)

        IF( coltab_dpv == -1 ) CALL xstctfn(coltabfn_dpv)
        CALL xsetclrs(coltab_dpv)
        CALL xcolor(1)    ! Set color to black
        CALL xcfont(3)
!
        call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
        cl(1)=0.0
        cl(2)=2.0           ! Set tentative contour interval
        CALL xnctrs(5,20)   ! Set lower and upper limits
        mode=2

        call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
        call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
        CALL xctrclr(ibgncol_dpv,iendcol_dpv)

        CALL xctrlim(0.0,30.0)

        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
        WRITE(6,'(a,f10.2)') &
       ' Plotting colfil Kdp field at elevation angle = ', eleva
        CALL xcolfil(vartilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

        IF ( overlayrefcontour == 1 ) THEN
          call XCLTYP(0)
          CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF

        CALL xwdwof
        CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

        CALL xchsiz(0.025*2*h_range_max )  ! Set character size
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

        WRITE(string,'(a,f5.2)') 'Kdp at Elevation Angle=',eleva
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
             h_range_max_n+2*h_range_max*0.04, trim(string) )

        WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcpalet(1)            ! Plot horizontal color palette

        xscale = h_range_max
        yscale = h_range_max
        xorig = -(h_range_max_e+h_range_max_w)/2.0
        yorig = -(h_range_max_s+h_range_max_n)/2.0

        xradar = 0.0
        yradar = 0.0
        CALL xthick(3)
        CALL xpenup( xradar-xscale*0.02, yradar )
        CALL xpendn( xradar+xscale*0.02, yradar )
        CALL xpenup( xradar, yradar-yscale*0.02 )
        CALL xpendn( xradar, yradar+yscale*0.02 )
        call xthick(1)
        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

        xscale = h_range_max_e+h_range_max_w
        yscale = h_range_max_s+h_range_max_n
        CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                    xscale,yscale,xorig,yorig)
        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        IF( ovrmap /= 0 ) THEN
          DO i=1,MIN(mxmapfile,nmapfile)
            lmapfile=len_trim(mapfile(i))
            WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

            INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
            IF( .NOT.fexist) THEN
              WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
            ELSE
              CALL xcolor(mapcol(i))
              IF(mapline_style(i) == 1) THEN
                CALL xthick(1)
                CALL xbrokn(6,3,6,3)
              ELSE IF(mapline_style(i) == 2) THEN
                CALL xthick(1)
              ELSE IF(mapline_style(i) == 3) THEN
                CALL xthick(3)
                CALL xfull
              END IF
              CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
            END IF
          END DO
        END IF
    
        CALL xwdwof
        CALL xframe
      END IF  ! pltkdp

      IF(wrtdualp > 0 .AND. pltrhohv > 0) THEN
        last_azim=-999.
        jazim=0
        refmax=-999.
        refmin=999.
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              vartilt(i,jazim)=max(rhohvvol(i,jidx,kref),0.)

              IF(rhohvvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(rhohvvol(i,jidx,kref),refmax)
                refmin=min(rhohvvol(i,jidx,kref),refmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        WRITE(6,'(a,2f10.2)') ' Rhohvmin, Rhohvmax on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
        CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        CALL xlabmask(1)

        IF( coltab_ref == -1 ) CALL xstctfn(coltabfn_dpv)
        CALL xsetclrs(coltab_dpv)
        CALL xcolor(1)    ! Set color to black
        CALL xcfont(3)
!
        call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
        cl(1)=0.0
        cl(2)=0.02           ! Set tentative contour interval
        CALL xnctrs(5,50)   ! Set lower and upper limits
        mode=2

        call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
        call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
        CALL xctrclr(ibgncol_ref,iendcol_ref)   ! Use colors between 10 and 39 inclusive

        CALL xctrlim(0.8,1.0) 

        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
        WRITE(6,'(a,f10.2)') &
       ' Plotting colfil rhohv field at elevation angle = ', eleva
        CALL xcolfil(reftilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

        IF ( overlayrefcontour == 1 ) THEN
          call XCLTYP(0)
          CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF

        CALL xwdwof
        CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

        CALL xchsiz(0.025*2*h_range_max )  ! Set character size
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

        WRITE(string,'(a,f5.2)') 'RhoHV at Elevation Angle=',eleva
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
             h_range_max_n+2*h_range_max*0.04, trim(string) )

        WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcpalet(1)            ! Plot horizontal color palette

        xscale = h_range_max
        yscale = h_range_max
        xorig = -(h_range_max_e+h_range_max_w)/2.0
        yorig = -(h_range_max_s+h_range_max_n)/2.0

        xradar = 0.0
        yradar = 0.0
        CALL xthick(3)
        CALL xpenup( xradar-xscale*0.02, yradar )
        CALL xpendn( xradar+xscale*0.02, yradar )
        CALL xpenup( xradar, yradar-yscale*0.02 )
        CALL xpendn( xradar, yradar+yscale*0.02 )
        call xthick(1)
        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

        xscale = h_range_max_e+h_range_max_w
        yscale = h_range_max_s+h_range_max_n
        CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                    xscale,yscale,xorig,yorig)
        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        IF( ovrmap /= 0 ) THEN
          DO i=1,MIN(mxmapfile,nmapfile)
            lmapfile=len_trim(mapfile(i))
            WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

            INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
            IF( .NOT.fexist) THEN
              WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
            ELSE
              CALL xcolor(mapcol(i))
              IF(mapline_style(i) == 1) THEN
                CALL xthick(1)
                CALL xbrokn(6,3,6,3)
              ELSE IF(mapline_style(i) == 2) THEN
                CALL xthick(1)
              ELSE IF(mapline_style(i) == 3) THEN
                CALL xthick(3)
                CALL xfull
              END IF
              CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
            END IF
          END DO
        END IF
    
        CALL xwdwof
        CALL xframe
      END IF  ! pltrhohv

      IF(wrtqx > 0 .AND. pltqr > 0) THEN
        last_azim=-999.
        jazim=0
        refmin = 1000.0
        refmax =-1000.0
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              vartilt(i,jazim)=min(max(qrvol(i,jidx,kref),0.),20.)

              IF(qrvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(qrvol(i,jidx,kref),refmax)
                refmin=min(qrvol(i,jidx,kref),refmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        WRITE(6,'(a,2f10.2)') ' qr min, qr max on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
        CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        CALL xlabmask(1)

        IF( coltab_dpv == -1 ) CALL xstctfn(coltabfn_dpv)
        CALL xsetclrs(coltab_dpv)
        CALL xcolor(1)    ! Set color to black
        CALL xcfont(3)
!
        call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
        cl(1)=0.0
        cl(2)=0.5           ! Set tentative contour interval
        CALL xnctrs(2,30)   ! Set lower and upper limits
        mode=2

        call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
        call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
        CALL xctrclr(ibgncol_dpv,iendcol_dpv)

        CALL xctrlim(-2.,5.)

        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
        WRITE(6,'(a,f10.2)') &
       ' Plotting colfil Zdr field at elevation angle = ', eleva
        CALL xcolfil(vartilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

        IF ( overlayrefcontour == 1 ) THEN
          call XCLTYP(0)
          CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF

        CALL xwdwof
        CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

        CALL xchsiz(0.025*2*h_range_max )  ! Set character size
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

        WRITE(string,'(a,f5.2)') 'Zdr at Elevation Angle=',eleva
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
             h_range_max_n+2*h_range_max*0.04, trim(string) )

        WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcpalet(1)            ! Plot horizontal color palette

        xscale = h_range_max
        yscale = h_range_max
        xorig = -(h_range_max_e+h_range_max_w)/2.0
        yorig = -(h_range_max_s+h_range_max_n)/2.0

        xradar = 0.0
        yradar = 0.0
        CALL xthick(3)
        CALL xpenup( xradar-xscale*0.02, yradar )
        CALL xpendn( xradar+xscale*0.02, yradar )
        CALL xpenup( xradar, yradar-yscale*0.02 )
        CALL xpendn( xradar, yradar+yscale*0.02 )
        call xthick(1)
        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

        xscale = h_range_max_e+h_range_max_w
        yscale = h_range_max_s+h_range_max_n
        CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                    xscale,yscale,xorig,yorig)
        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        IF( ovrmap /= 0 ) THEN
          DO i=1,MIN(mxmapfile,nmapfile)
            lmapfile=len_trim(mapfile(i))
            WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

            INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
            IF( .NOT.fexist) THEN
              WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
            ELSE
              CALL xcolor(mapcol(i))
              IF(mapline_style(i) == 1) THEN
                CALL xthick(1)
                CALL xbrokn(6,3,6,3)
              ELSE IF(mapline_style(i) == 2) THEN
                CALL xthick(1)
              ELSE IF(mapline_style(i) == 3) THEN
                CALL xthick(3)
                CALL xfull
              END IF
              CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
            END IF
          END DO
        END IF
    
        CALL xwdwof
        CALL xframe
      END IF  ! pltqr

      IF(wrtqx > 0 .AND. pltqs > 0) THEN
        last_azim=-999.
        jazim=0
        refmin = 1000.0
        refmax =-1000.0
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              vartilt(i,jazim)=min(max(qsvol(i,jidx,kref),0.),20.)

              IF(qsvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(qsvol(i,jidx,kref),refmax)
                refmin=min(qsvol(i,jidx,kref),refmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        WRITE(6,'(a,2f10.2)') ' qs min, qs max on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
        CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        CALL xlabmask(1)

        IF( coltab_dpv == -1 ) CALL xstctfn(coltabfn_dpv)
        CALL xsetclrs(coltab_dpv)
        CALL xcolor(1)    ! Set color to black
        CALL xcfont(3)
!
        call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
        cl(1)=0.0
        cl(2)=0.5           ! Set tentative contour interval
        CALL xnctrs(2,30)   ! Set lower and upper limits
        mode=2

        call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
        call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
        CALL xctrclr(ibgncol_dpv,iendcol_dpv)

        CALL xctrlim(-2.,5.)

        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
        WRITE(6,'(a,f10.2)') &
       ' Plotting colfil Zdr field at elevation angle = ', eleva
        CALL xcolfil(vartilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

        IF ( overlayrefcontour == 1 ) THEN
          call XCLTYP(0)
          CALL xctrclr(14,15)
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF

        CALL xwdwof
        CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

        CALL xchsiz(0.025*2*h_range_max )  ! Set character size
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

        WRITE(string,'(a,f5.2)') 'Zdr at Elevation Angle=',eleva
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
             h_range_max_n+2*h_range_max*0.04, trim(string) )

        WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcpalet(1)            ! Plot horizontal color palette

        xscale = h_range_max
        yscale = h_range_max
        xorig = -(h_range_max_e+h_range_max_w)/2.0
        yorig = -(h_range_max_s+h_range_max_n)/2.0

        xradar = 0.0
        yradar = 0.0
        CALL xthick(3)
        CALL xpenup( xradar-xscale*0.02, yradar )
        CALL xpendn( xradar+xscale*0.02, yradar )
        CALL xpenup( xradar, yradar-yscale*0.02 )
        CALL xpendn( xradar, yradar+yscale*0.02 )
        call xthick(1)
        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

        xscale = h_range_max_e+h_range_max_w
        yscale = h_range_max_s+h_range_max_n
        CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                    xscale,yscale,xorig,yorig)
        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        IF( ovrmap /= 0 ) THEN
          DO i=1,MIN(mxmapfile,nmapfile)
            lmapfile=len_trim(mapfile(i))
            WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

            INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
            IF( .NOT.fexist) THEN
              WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
            ELSE
              CALL xcolor(mapcol(i))
              IF(mapline_style(i) == 1) THEN
                CALL xthick(1)
                CALL xbrokn(6,3,6,3)
              ELSE IF(mapline_style(i) == 2) THEN
                CALL xthick(1)
              ELSE IF(mapline_style(i) == 3) THEN
                CALL xthick(3)
                CALL xfull
              END IF
              CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
            END IF
          END DO
        END IF
    
        CALL xwdwof
        CALL xframe
      END IF  ! pltqs

      IF(wrtdualp > 0 .AND. pltqh > 0) THEN
        last_azim=-999.
        jazim=0
        refmin = 1000.0
        refmax =-1000.0
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              vartilt(i,jazim)=min(max(qhvol(i,jidx,kref),0.),20.)

              IF(qhvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(qhvol(i,jidx,kref),refmax)
                refmin=min(qhvol(i,jidx,kref),refmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        WRITE(6,'(a,2f10.2)') ' qh min, qh max on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
        CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        CALL xlabmask(1)

        IF( coltab_dpv == -1 ) CALL xstctfn(coltabfn_dpv)
        CALL xsetclrs(coltab_dpv)
        CALL xcolor(1)    ! Set color to black
        CALL xcfont(3)
!
        call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
        cl(1)=0.0
        cl(2)=0.5           ! Set tentative contour interval
        CALL xnctrs(2,30)   ! Set lower and upper limits
        mode=2

        call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
        call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
        CALL xctrclr(ibgncol_dpv,iendcol_dpv)

        CALL xctrlim(-2.,5.)

        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
        WRITE(6,'(a,f10.2)') &
       ' Plotting colfil Zdr field at elevation angle = ', eleva
        CALL xcolfil(vartilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

        IF ( overlayrefcontour == 1 ) THEN
          call XCLTYP(0)
          CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF

        CALL xwdwof
        CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

        CALL xchsiz(0.025*2*h_range_max )  ! Set character size
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

        WRITE(string,'(a,f5.2)') 'Zdr at Elevation Angle=',eleva
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
             h_range_max_n+2*h_range_max*0.04, trim(string) )

        WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcpalet(1)            ! Plot horizontal color palette

        xscale = h_range_max
        yscale = h_range_max
        xorig = -(h_range_max_e+h_range_max_w)/2.0
        yorig = -(h_range_max_s+h_range_max_n)/2.0

        xradar = 0.0
        yradar = 0.0
        CALL xthick(3)
        CALL xpenup( xradar-xscale*0.02, yradar )
        CALL xpendn( xradar+xscale*0.02, yradar )
        CALL xpenup( xradar, yradar-yscale*0.02 )
        CALL xpendn( xradar, yradar+yscale*0.02 )
        call xthick(1)
        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

        xscale = h_range_max_e+h_range_max_w
        yscale = h_range_max_s+h_range_max_n
        CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                    xscale,yscale,xorig,yorig)
        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        IF( ovrmap /= 0 ) THEN
          DO i=1,MIN(mxmapfile,nmapfile)
            lmapfile=len_trim(mapfile(i))
            WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

            INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
            IF( .NOT.fexist) THEN
              WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
            ELSE
              CALL xcolor(mapcol(i))
              IF(mapline_style(i) == 1) THEN
                CALL xthick(1)
                CALL xbrokn(6,3,6,3)
              ELSE IF(mapline_style(i) == 2) THEN
                CALL xthick(1)
              ELSE IF(mapline_style(i) == 3) THEN
                CALL xthick(3)
                CALL xfull
              END IF
              CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
            END IF
          END DO
        END IF
    
        CALL xwdwof
        CALL xframe
      END IF  ! pltqh

      IF(wrtuvwt > 0 .AND. pltu > 0) THEN
 
        cl(1)=0.0
        cl(2)=0.5           ! Set tentative contour interval

        last_azim=-999.
        jazim=0
        refmin = 1000.0
        refmax =-1000.0
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              vartilt(i,jazim)=min(max(uvol(i,jidx,kref),-90.),90.)

              IF(uvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(uvol(i,jidx,kref),refmax)
                refmin=min(uvol(i,jidx,kref),refmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        WRITE(6,'(a,2f10.2)') ' u min, u max on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
        CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        CALL xlabmask(1)

        IF( coltab_dpv == -1 ) CALL xstctfn(coltabfn_dpv)
        CALL xsetclrs(coltab_dpv)
        CALL xcolor(1)    ! Set color to black
        CALL xcfont(3)
!
        call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
        cl(1)=0.0
        cl(2)=0.5           ! Set tentative contour interval
        CALL xnctrs(2,30)   ! Set lower and upper limits
        mode=2

        call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
        call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
        CALL xctrclr(ibgncol_dpv,iendcol_dpv)

        CALL xctrlim(-2.,5.)

        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
        WRITE(6,'(a,f10.2)') &
       ' Plotting colfil Zdr field at elevation angle = ', eleva
        CALL xcolfil(vartilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

        IF ( overlayrefcontour == 1 ) THEN
          call XCLTYP(0)
          CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF

        CALL xwdwof
        CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

        CALL xchsiz(0.025*2*h_range_max )  ! Set character size
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

        WRITE(string,'(a,f5.2)') 'Zdr at Elevation Angle=',eleva
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
             h_range_max_n+2*h_range_max*0.04, trim(string) )

        WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcpalet(1)            ! Plot horizontal color palette

        xscale = h_range_max
        yscale = h_range_max
        xorig = -(h_range_max_e+h_range_max_w)/2.0
        yorig = -(h_range_max_s+h_range_max_n)/2.0

        xradar = 0.0
        yradar = 0.0
        CALL xthick(3)
        CALL xpenup( xradar-xscale*0.02, yradar )
        CALL xpendn( xradar+xscale*0.02, yradar )
        CALL xpenup( xradar, yradar-yscale*0.02 )
        CALL xpendn( xradar, yradar+yscale*0.02 )
        call xthick(1)
        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

        xscale = h_range_max_e+h_range_max_w
        yscale = h_range_max_s+h_range_max_n
        CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                    xscale,yscale,xorig,yorig)
        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        IF( ovrmap /= 0 ) THEN
          DO i=1,MIN(mxmapfile,nmapfile)
            lmapfile=len_trim(mapfile(i))
            WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

            INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
            IF( .NOT.fexist) THEN
              WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
            ELSE
              CALL xcolor(mapcol(i))
              IF(mapline_style(i) == 1) THEN
                CALL xthick(1)
                CALL xbrokn(6,3,6,3)
              ELSE IF(mapline_style(i) == 2) THEN
                CALL xthick(1)
              ELSE IF(mapline_style(i) == 3) THEN
                CALL xthick(3)
                CALL xfull
              END IF
              CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
            END IF
          END DO
        END IF
    
        CALL xwdwof
        CALL xframe
      END IF  ! pltu

      IF(wrtuvwt > 0 .AND. pltv > 0) THEN
 
        cl(1)=0.0
        cl(2)=1.0           ! Set tentative contour interval

        last_azim=-999.
        jazim=0
        refmin = 1000.0
        refmax =-1000.0
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              vartilt(i,jazim)=min(max(vvol(i,jidx,kref),-90.),90.)

              IF(vvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(vvol(i,jidx,kref),refmax)
                refmin=min(vvol(i,jidx,kref),refmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        WRITE(6,'(a,2f10.2)') '  v min,  v max on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
        CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        CALL xlabmask(1)

        IF( coltab_dpv == -1 ) CALL xstctfn(coltabfn_dpv)
        CALL xsetclrs(coltab_dpv)
        CALL xcolor(1)    ! Set color to black
        CALL xcfont(3)
!
        call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
        cl(1)=0.0
        cl(2)=0.5           ! Set tentative contour interval
        CALL xnctrs(2,30)   ! Set lower and upper limits
        mode=2

        call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
        call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
        CALL xctrclr(ibgncol_dpv,iendcol_dpv)

        CALL xctrlim(-2.,5.)

        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
        WRITE(6,'(a,f10.2)') &
       ' Plotting colfil Zdr field at elevation angle = ', eleva
        CALL xcolfil(vartilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

        IF ( overlayrefcontour == 1 ) THEN
          call XCLTYP(0)
          CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF

        CALL xwdwof
        CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

        CALL xchsiz(0.025*2*h_range_max )  ! Set character size
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

        WRITE(string,'(a,f5.2)') 'Zdr at Elevation Angle=',eleva
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
             h_range_max_n+2*h_range_max*0.04, trim(string) )

        WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcpalet(1)            ! Plot horizontal color palette

        xscale = h_range_max
        yscale = h_range_max
        xorig = -(h_range_max_e+h_range_max_w)/2.0
        yorig = -(h_range_max_s+h_range_max_n)/2.0

        xradar = 0.0
        yradar = 0.0
        CALL xthick(3)
        CALL xpenup( xradar-xscale*0.02, yradar )
        CALL xpendn( xradar+xscale*0.02, yradar )
        CALL xpenup( xradar, yradar-yscale*0.02 )
        CALL xpendn( xradar, yradar+yscale*0.02 )
        call xthick(1)
        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

        xscale = h_range_max_e+h_range_max_w
        yscale = h_range_max_s+h_range_max_n
        CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                    xscale,yscale,xorig,yorig)
        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        IF( ovrmap /= 0 ) THEN
          DO i=1,MIN(mxmapfile,nmapfile)
            lmapfile=len_trim(mapfile(i))
            WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

            INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
            IF( .NOT.fexist) THEN
              WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
            ELSE
              CALL xcolor(mapcol(i))
              IF(mapline_style(i) == 1) THEN
                CALL xthick(1)
                CALL xbrokn(6,3,6,3)
              ELSE IF(mapline_style(i) == 2) THEN
                CALL xthick(1)
              ELSE IF(mapline_style(i) == 3) THEN
                CALL xthick(3)
                CALL xfull
              END IF
              CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
            END IF
          END DO
        END IF
    
        CALL xwdwof
        CALL xframe
      END IF  ! pltv

      IF(wrtuvwt > 0 .AND. pltw > 0) THEN
        last_azim=-999.
        jazim=0
        refmin = 1000.0
        refmax =-1000.0
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              vartilt(i,jazim)=min(max(wvol(i,jidx,kref),-90.),90.)

              IF(wvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(wvol(i,jidx,kref),refmax)
                refmin=min(wvol(i,jidx,kref),refmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        WRITE(6,'(a,2f10.2)') ' Zdrmin, Zdrmax on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
        CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        CALL xlabmask(1)

        IF( coltab_dpv == -1 ) CALL xstctfn(coltabfn_dpv)
        CALL xsetclrs(coltab_dpv)
        CALL xcolor(1)    ! Set color to black
        CALL xcfont(3)
!
        call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
        cl(1)=0.0
        cl(2)=0.5           ! Set tentative contour interval
        CALL xnctrs(2,30)   ! Set lower and upper limits
        mode=2

        call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
        call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
        CALL xctrclr(ibgncol_dpv,iendcol_dpv)

        CALL xctrlim(-2.,5.)

        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
        WRITE(6,'(a,f10.2)') &
       ' Plotting colfil Zdr field at elevation angle = ', eleva
        CALL xcolfil(vartilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

        IF ( overlayrefcontour == 1 ) THEN
          call XCLTYP(0)
          CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF

        CALL xwdwof
        CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

        CALL xchsiz(0.025*2*h_range_max )  ! Set character size
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

        WRITE(string,'(a,f5.2)') 'Zdr at Elevation Angle=',eleva
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
             h_range_max_n+2*h_range_max*0.04, trim(string) )

        WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcpalet(1)            ! Plot horizontal color palette

        xscale = h_range_max
        yscale = h_range_max
        xorig = -(h_range_max_e+h_range_max_w)/2.0
        yorig = -(h_range_max_s+h_range_max_n)/2.0

        xradar = 0.0
        yradar = 0.0
        CALL xthick(3)
        CALL xpenup( xradar-xscale*0.02, yradar )
        CALL xpendn( xradar+xscale*0.02, yradar )
        CALL xpenup( xradar, yradar-yscale*0.02 )
        CALL xpendn( xradar, yradar+yscale*0.02 )
        call xthick(1)
        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

        xscale = h_range_max_e+h_range_max_w
        yscale = h_range_max_s+h_range_max_n
        CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                    xscale,yscale,xorig,yorig)
        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        IF( ovrmap /= 0 ) THEN
          DO i=1,MIN(mxmapfile,nmapfile)
            lmapfile=len_trim(mapfile(i))
            WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

            INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
            IF( .NOT.fexist) THEN
              WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
            ELSE
              CALL xcolor(mapcol(i))
              IF(mapline_style(i) == 1) THEN
                CALL xthick(1)
                CALL xbrokn(6,3,6,3)
              ELSE IF(mapline_style(i) == 2) THEN
                CALL xthick(1)
              ELSE IF(mapline_style(i) == 3) THEN
                CALL xthick(3)
                CALL xfull
              END IF
              CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
            END IF
          END DO
        END IF
    
        CALL xwdwof
        CALL xframe
      END IF  ! pltw

      IF(wrtuvwt > 0 .AND. pltt > 0) THEN
        last_azim=-999.
        jazim=0
        refmin = 1000.0
        refmax =-1000.0
        DO j=1,maxazim
          jidx=j
          IF(jidx > maxazim) jidx=jidx-maxazim
          data_azim=azmrvol(jidx,kref)
          IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
            jazim=jazim+1
            last_azim=data_azim
            xconst=cos( (90.-data_azim)*deg2rad ) 
            yconst=sin( (90.-data_azim)*deg2rad ) 
            DO i=1,maxrefgate
              rad = 0.001 * rngrvol(i,kref)
              xrplt(i,jazim)= rad * xconst
              yrplt(i,jazim)= rad * yconst
              vartilt(i,jazim)=min(max(tvol(i,jidx,kref),-40.),50.)

              IF(tvol(i,jidx,kref) > missing_value ) THEN 
                refmax=max(tvol(i,jidx,kref),refmax)
                refmin=min(tvol(i,jidx,kref),refmin)
              END IF
            END DO
          END IF
        END DO
        nazim=jazim

        refmin=refmin-273.15
        refmax=refmax-273.15
        WRITE(6,'(a,2f10.2)') ' T min, T max on tilt= ', refmin, refmax
        WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
        CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
        CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        CALL xlabmask(1)

        IF( coltab_dpv == -1 ) CALL xstctfn(coltabfn_dpv)
        CALL xsetclrs(coltab_dpv)
        CALL xcolor(1)    ! Set color to black
        CALL xcfont(3)
!
        call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
        cl(1)=0.0
        cl(2)=0.5           ! Set tentative contour interval
        CALL xnctrs(2,30)   ! Set lower and upper limits
        mode=2

        call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
        call xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
        CALL xctrclr(ibgncol_dpv,iendcol_dpv)

        CALL xctrlim(-2.,5.)

        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
        WRITE(6,'(a,f10.2)') &
       ' Plotting colfil T field at elevation angle = ', eleva
        CALL xcolfil(vartilt,xrplt,yrplt,irwxy,xrw,yrw,                   &
                       maxrefgate,maxrefgate,nazim, cl, ncl,mode)

        IF ( overlayrefcontour == 1 ) THEN
          call XCLTYP(0)
          CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy, &
                    maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF

        CALL xwdwof
        CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

        CALL xchsiz(0.025*2*h_range_max )  ! Set character size
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

        WRITE(string,'(a,f5.2)') 'Zdr at Elevation Angle=',eleva
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
             h_range_max_n+2*h_range_max*0.04, trim(string) )

        WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
        CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcpalet(1)            ! Plot horizontal color palette

        xscale = h_range_max
        yscale = h_range_max
        xorig = -(h_range_max_e+h_range_max_w)/2.0
        yorig = -(h_range_max_s+h_range_max_n)/2.0

        xradar = 0.0
        yradar = 0.0
        CALL xthick(3)
        CALL xpenup( xradar-xscale*0.02, yradar )
        CALL xpendn( xradar+xscale*0.02, yradar )
        CALL xpenup( xradar, yradar-yscale*0.02 )
        CALL xpendn( xradar, yradar+yscale*0.02 )
        call xthick(1)
        CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
        CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

        xscale = h_range_max_e+h_range_max_w
        yscale = h_range_max_s+h_range_max_n
        CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon, &
                    xscale,yscale,xorig,yorig)
        CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

        IF( ovrmap /= 0 ) THEN
          DO i=1,MIN(mxmapfile,nmapfile)
            lmapfile=len_trim(mapfile(i))
            WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

            INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
            IF( .NOT.fexist) THEN
              WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
            ELSE
              CALL xcolor(mapcol(i))
              IF(mapline_style(i) == 1) THEN
                CALL xthick(1)
                CALL xbrokn(6,3,6,3)
              ELSE IF(mapline_style(i) == 2) THEN
                CALL xthick(1)
              ELSE IF(mapline_style(i) == 3) THEN
                CALL xthick(3)
                CALL xfull
              END IF
              CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
            END IF
          END DO
        END IF
    
        CALL xwdwof
        CALL xframe
      END IF  ! pltt

    END DO ! ktilt 

    DEALLOCATE( xrplt )
    DEALLOCATE( yrplt )
    DEALLOCATE( reftilt )
    DEALLOCATE( vartilt )
    DEALLOCATE( irwxy )
    DEALLOCATE( xrw )
    DEALLOCATE( yrw )

  END IF  !  n_elevation_ref>0
!
!-----------------------------------------------------------------------
! To plot vertical cross-sections of reflectivity
!-----------------------------------------------------------------------
!
  IF(rhimode > 0 .OR. n_azimuth_ref > 0) THEN

    print *, ' maxrefgate= ',maxrefgate
    print *, ' nref_tilts= ',nref_tilts
    print *, ' n_azimuth_ref = ',n_azimuth_ref
    ALLOCATE( rplt(maxrefgate,nref_tilts) , stat=ialstatus)
    ALLOCATE( rpltm(maxrefgate,nref_tilts), stat=ialstatus )
    ALLOCATE( zplt(maxrefgate,nref_tilts), stat=ialstatus )
    ALLOCATE( datarz(maxrefgate,nref_tilts), stat=ialstatus )
    ALLOCATE( iwrz(maxrefgate,nref_tilts), stat=ialstatus )
    ALLOCATE( xrw(8*maxrefgate), stat=ialstatus )
    ALLOCATE( yrw(8*maxrefgate), stat=ialstatus )

    zradar=0.001*radar_elevation
    DO k=1,nref_tilts
      DO i=1,maxrefgate
        CALL beamhgt(elvref(k),rngrvol(i,k),height,sfcrng)
        rplt(i,k)=0.001*sfcrng
        rpltm(i,k)=-0.001*sfcrng
        zplt(i,k)=zradar+0.001*height
      END DO
    END DO

    xtick=0.1*h_range_max_xs
    xmax=h_range_max_xs

    ztick=0.5
    zmin=0.0
    zmax=z_range_max

    IF(rhidata > 0) THEN
      n_azimuth_ref=1
      azimuth_ref(1)=azmrvol(1,1)
    END IF

    IF(rhimode > 0) THEN
      xmin=0.0
      xlab=0.5*h_range_max_xs
    ELSE
      xmin=-h_range_max_xs
      xlab=0.0
    END IF

    DO iazimuth= 1,n_azimuth_ref

      WRITE(6,'(a,f10.2)') &
        ' Plotting vertical cross-sections at azimuth angle = ',      &
                                azimuth_ref(iazimuth)

      IF(pltref > 0) THEN

        varname='Reflectivity (dBZ)'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=5.0

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,refvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_ref,ibgncol_ref,iendcol_ref,coltabfn_ref,      &
              pltmin_ref,pltmax_ref,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

      END IF ! pltref

      IF(wrtdualp > 0 .AND. pltzdr > 0) THEN

        varname='Zdr (dB)'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=0.5

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,zdrvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              -2.0,6.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

        varname='Hail Flag'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=0.05

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,hailflg,&
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              0.0,1.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

      END IF ! pltzdr

      IF(wrtdualp > 0 .AND. pltkdp > 0) THEN

        varname='Kdp (degrees/km)'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=2.0

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,kdpvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              0.0,30.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

      END IF ! pltkdp

      IF(wrtdualp > 0 .AND. pltrhohv > 0) THEN

        varname='rhoHV'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=0.05

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,rhohvvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              0.0,1.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

      END IF ! rhohvplt

      IF(wrtqx > 0 .AND. pltqr > 0) THEN
 
        varname='Rain Mixing Ratio (g/kg)'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=1.0

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,qrvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              0.0,20.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

      END IF ! pltqr

      IF(wrtqx > 0 .AND. pltqs > 0) THEN
 
        varname='Snow Mixing Ratio (g/kg)'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=0.1

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,qsvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              0.0,2.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

      END IF ! pltqs

      IF(wrtqx > 0 .AND. pltqh > 0) THEN
 
        varname='Hail Mixing Ratio (g/kg)'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=1.0

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,qhvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              0.0,20.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

      END IF ! pltqh

      IF(wrtuvwt > 0 .AND. pltu > 0) THEN
 
        varname='u Wind Component (m/s)'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=5.0

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,uvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              -50.0,50.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

      END IF ! uplt

      IF(wrtuvwt > 0 .AND. pltv > 0) THEN
 
        varname='v Wind Component (m/s)'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=5.0

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,vvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              -50.0,50.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)
 
      END IF ! vplt

      IF(wrtuvwt > 0 .AND. pltw > 0) THEN
 
        varname='Vertical Velocity (m/s)'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=0.5

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,wvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              -10.0,10.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

      END IF ! pltw

      IF(wrtuvwt > 0 .AND. pltt > 0) THEN
 
        varname='Air Temperature (C)'
        WRITE(6,'(a,a)') ' Plotting ',TRIM(varname)
        cl=0.0
        cl(2)=5.0

        CALL pltradxs(maxrefgate,maxazim,nref_tilts,azmrvol,tvol, &
              rplt,rpltm,zplt,datarz,                               &
              varname,azimuth_ref(iazimuth),rhimode,missing_value,  &
              xmin,xmax,xtick,zmin,zmax,ztick,                      &
              coltab_dpv,ibgncol_dpv,iendcol_dpv,coltabfn_dpv,      &
              -65.0,40.0,iyear,imon,iday,ihour,imin,     &
              cl,iwrz,xrw,yrw)

      END IF ! pltt

    END DO ! iazimuth

    DEALLOCATE( rplt )
    DEALLOCATE( rpltm )
    DEALLOCATE( zplt )
    DEALLOCATE( datarz )
    DEALLOCATE( iwrz )
    DEALLOCATE( xrw )
    DEALLOCATE( yrw )

  END IF   ! n_azimuth_ref
!
! Velocity data plotting section
!
  WRITE(6,'(a)') ' Velocity Plotting Section'
  IF (n_elevation_vel > 0 ) THEN

    WRITE(6,'(a)') ' Allocating x,y,veltilt ...'
    ALLOCATE( xvplt(maxvelgate,(maxazim+1)), STAT=ialloc)
    ALLOCATE( yvplt(maxvelgate,(maxazim+1)), STAT=ialloc )
    ALLOCATE( veltilt(maxvelgate,(maxazim+1)), STAT=ialloc )
    ALLOCATE( ivwxy(maxvelgate,(maxazim+1)), STAT=ialloc )
    ALLOCATE( xvw(8*maxvelgate), STAT=ialloc )
    ALLOCATE( yvw(8*maxvelgate), STAT=ialloc )

    IF ( overlayrefcontour == 1 ) THEN
      ALLOCATE( xrplt(maxrefgate,(maxazim+1)), STAT=ialloc )
      ALLOCATE( yrplt(maxrefgate,(maxazim+1)), STAT=ialloc )
      ALLOCATE( reftilt(maxrefgate,(maxazim+1)), STAT=ialloc )
      ALLOCATE( irwxy(maxrefgate,(maxazim+1)), STAT=ialloc )
      ALLOCATE( xrw(8*maxrefgate), STAT=ialloc )
      ALLOCATE( yrw(8*maxrefgate), STAT=ialloc )
    END IF

    WRITE(6,'(a)') ' Radial Velocity tilt plotting...'

    DO ktilt = 1,n_elevation_vel

      elev_to_plot = elevations_vel(ktilt)

      kvel = 0 
      DO k = 1,maxelev
        IF( abs( elvref(k)-elev_to_plot )                   &
                            < elev_error_range ) THEN
          kvel = k
          EXIT
        END IF
      END DO

      IF(kvel == 0 ) CYCLE  ! no elevation angle 
                            ! within error range was found
      velmin = 1000.0
      velmax =-1000.0
      velmin = 1000.0
      velmax =-1000.0

      eleva=elvref(kvel)

      last_azim=-999.
      jazim=0
      DO j=1,maxazim+1
        jidx=j
        IF(jidx > maxazim) jidx=jidx-maxazim
        data_azim=azmvvol(jidx,kvel)
        IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
          jazim=jazim+1
          last_azim=data_azim
          xconst=cos( (90.-data_azim)*deg2rad ) 
          yconst=sin( (90.-data_azim)*deg2rad ) 
          DO i=1,maxvelgate
            rad = 0.001 * rngvvol(i,kvel)
            xvplt(i,jazim)= rad * xconst
            yvplt(i,jazim)= rad * yconst
            veltilt(i,jazim)=max(velvol(i,jidx,kvel),-150.)

            IF(velvol(i,jidx,kvel) > missing_value ) THEN 
              velmax=max(velvol(i,jidx,kvel),velmax)
              velmin=min(velvol(i,jidx,kvel),velmin)
            END IF
          END DO
        END IF
      END DO

      nazim=jazim

      WRITE(6,'(a,2f10.2)') ' Velmin, Velmax on tilt= ', velmin, velmax
      WRITE(6,'(a,i5)') 'N unique azimuths: ',nazim
!
!-----------------------------------------------------------------------
! Define plotting space.
!
      CALL xpspac(0.5-factorx,0.5+factorx,0.5-factory,0.5+factory)
      CALL xmap(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

      CALL xlabmask(1)

      IF( coltab_vel == -1 ) CALL xstctfn(coltabfn_vel)
      CALL xsetclrs(coltab_vel)  ! Use internally defined color map number 6
      CALL xcolor(1)    ! Set color to black
      CALL xcfont(3)
!
      call xaxfmt('(I4)')
!
! Set tentative color interval and limit on the min/max number of contours
!
      cl(1)=0.0
      cl(2)=5.0           ! Set tentative contour interval
      CALL xnctrs(2,30)  ! Set lower and upper limits
                          ! on the number of contours allowed
                          ! It's used by both xcolfil and xclimit
      mode=2

      call xctrbadv(1)    ! Turn on bad value checking
                          !     in contouring routine
      call xbadval( missing_value ) ! Specify bad value flag

      CALL xctrclr(ibgncol_vel,iendcol_vel)

      CALL xctrlim(pltmin_vel,pltmax_vel)

      CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)
      WRITE(6,'(a,f10.2)') &
        'Plotting colfil radial velocity at elevation angle=',eleva
      CALL xcolfil(veltilt,xvplt,yvplt,ivwxy,xvw,yvw,                   &
                       maxvelgate,maxvelgate,nazim, cl, ncl,mode)

      IF ( overlayrefcontour == 1 ) THEN

        kref = 0 
        DO k = 1,maxelev
          IF( abs( elvref(k)-elev_to_plot )                   &
                              < elev_error_range ) THEN
            kref = k
            EXIT
          END IF
        END DO

        IF(kref /= 0 ) THEN
          refmin = 1000.0
          refmax =-1000.0
          velmin = 1000.0
          refmax =-1000.0

          eleva=elvref(kref)

          last_azim=-999.
          jazim=0
          DO j=1,maxazim+1
            jidx=j
            IF(jidx > maxazim) jidx=jidx-maxazim
            data_azim=azmrvol(jidx,kref)
            IF(data_azim >= 0. .AND. (abs(data_azim-last_azim) > 1.0E-04)) THEN
              jazim=jazim+1
              last_azim=data_azim
              xconst=cos( (90.-data_azim)*deg2rad ) 
              yconst=sin( (90.-data_azim)*deg2rad ) 
              DO i=1,maxrefgate
                rad = 0.001 * rngrvol(i,kref)
                xrplt(i,jazim)= rad * xconst
                yrplt(i,jazim)= rad * yconst
                reftilt(i,jazim)=max(refvol(i,jidx,kref),-30.)
    
                IF(refvol(i,jidx,kref) > missing_value ) THEN 
                  refmax=max(refvol(i,jidx,kref),refmax)
                  refmin=min(refvol(i,jidx,kref),refmin)
                END IF
              END DO
            END IF
          END DO
          nazim=jazim
    
          call XCLTYP(0)
          CALL xctrclr(14,15)   ! Use colors between 10 and 39 inclusive
          mode=3
          ncl=1
          cl(1)=15.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy,    &
                      maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=30.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy,    &
                      maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=45.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy,    &
                      maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
          cl(1)=55.0
          CALL XCONTA(reftilt,xrplt,yrplt,irwxy,    &
                      maxrefgate,maxrefgate,nazim, CL, NCL,MODE)
        END IF
      END IF

      CALL xwdwof
      CALL xaxsca(-h_range_max_w, h_range_max_e, 0.1*h_range_max,          &
                  -h_range_max_s, h_range_max_n, 0.1*h_range_max)

      CALL xchsiz(0.025*2*h_range_max )  ! Set character size
      CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

      WRITE(string,'(a,f5.2)') 'Radial Velocity at Elevation Angle=',eleva
      CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.04, trim(string) )

      WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'Scan Time:',      &
            imon,'/',iday,'/',iyear, ' ',ihour,':',imin 
      CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
           h_range_max_n+2*h_range_max*0.01, trim(stringtime) )

      CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
      CALL xcpalet(1)            ! Plot horizontal color palette

      xscale = h_range_max
      yscale = h_range_max
      xorig = -(h_range_max_e+h_range_max_w)/2.0
      yorig = -(h_range_max_s+h_range_max_n)/2.0

      xradar = 0.0
      yradar = 0.0
      CALL xthick(3)
      CALL xpenup( xradar-xscale*0.02, yradar )
      CALL xpendn( xradar+xscale*0.02, yradar )
      CALL xpenup( xradar, yradar-yscale*0.02 )
      CALL xpendn( xradar, yradar+yscale*0.02 )
      call xthick(1)
      CALL xchsiz(0.02* 2*h_range_max )  ! Set character size
      CALL xcharc(xradar, yradar-0.05*yscale, radar_symbol)

      xscale = h_range_max_e+h_range_max_w
      yscale = h_range_max_s+h_range_max_n
      CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,radar_lat,radar_lon,& 
                    xscale,yscale,xorig,yorig)
      CALL xwindw(-h_range_max_w, h_range_max_e, -h_range_max_s, h_range_max_n)

      IF( ovrmap /= 0 ) THEN
        DO i=1,MIN(mxmapfile,nmapfile)
          lmapfile=len_trim(mapfile(i))
          WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

          INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
          IF( .NOT.fexist) THEN
            WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                 ' not found. Corresponding map no plotted.'
          ELSE
            CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
          END IF
        END DO
      END IF

      CALL xwdwof

      CALL xframe

    END DO ! ktilt 

    DEALLOCATE( ivwxy )
    DEALLOCATE( xvw )
    DEALLOCATE( yvw )
    DEALLOCATE( xvplt )
    DEALLOCATE( yvplt )
    DEALLOCATE( veltilt )

    IF(overlayrefcontour == 1) THEN
      DEALLOCATE( irwxy )
      DEALLOCATE( xrw )
      DEALLOCATE( yrw )
      DEALLOCATE( xrplt )
      DEALLOCATE( yrplt )
      DEALLOCATE( reftilt )
    END IF

  END IF  !  n_elevation_vel>0
!
!-----------------------------------------------------------------------
! To plot vertical cross-sections of radial velocity
!-----------------------------------------------------------------------
!
  IF(n_azimuth_vel > 0) THEN

    ALLOCATE( rplt(maxvelgate,nvel_tilts))
    ALLOCATE( rpltm(maxvelgate,nvel_tilts))
    ALLOCATE( zplt(maxvelgate,nvel_tilts))
    ALLOCATE( datarz(maxvelgate,nvel_tilts) )
    ALLOCATE( iwrz(maxvelgate,nvel_tilts) )
    ALLOCATE( xvw(8*maxvelgate) )
    ALLOCATE( yvw(8*maxvelgate) )

    DO k=1,nvel_tilts
      DO i=1,maxvelgate
        rad = 0.001 * rngvvol(i,k)
        rplt(i,k)= rad*cos( elvref(k)*deg2rad )
        rpltm(i,k)=-rad*cos( elvref(k)*deg2rad )
        zplt(i,k)= rad*sin( elvref(k)*deg2rad )
      END DO
    END DO

    DO iazimuth= 1,n_azimuth_vel

      WRITE(6,'(a,f10.2)') &
        ' Plotting vertical cross-section at azimuth angle = ',      &
                                azimuth_vel(iazimuth)

      CALL xpspac(0.1,0.9,0.2,0.6)
      z_range_max = 12.0

      CALL xmap(-h_range_max, h_range_max, 0.0, z_range_max)

      CALL xlabmask(1)

      IF( coltab_vel == -1 ) CALL xstctfn(coltabfn_vel)
      CALL xsetclrs(coltab_vel)  ! Use internally defined color map number 6
      CALL xcolor(1)    ! Set color to black
      CALL xcfont(3)
!
! Set tentative color interval and limit on the min/max number of contours
!
      cl(1)=0.0
      cl(2)=5.0           ! Set tentative contour interval
      CALL xnctrs(2,100)  ! Set lower and upper limits
                          ! on the number of contours allowed
                          ! It's used by both xcolfil and xclimit
      CALL xnctrs(5,20)  ! Set lower and upper limits
      mode=2

      CALL xctrbadv(1)   ! Turn on bad value checking
                         ! in contouring routine
      CALL xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field for positive values
!
      CALL xctrclr(ibgncol_vel,iendcol_vel) ! Use colors between 10 and 39 inclusive

      CALL xctrlim(pltmin_vel,pltmax_vel) ! Only contours above zero are plotted
      CALL xwindw(-h_range_max, h_range_max, 0.0, z_range_max)

      DO k=1,nvel_tilts
        jj=1
        dazmin=720.
        DO j=1,maxazim
          delazm=azmrvol(j,k)-azimuth_vel(iazimuth)
          IF(delazm > 180.) delazm=delazm-360.
          IF(delazm < -180.) delazm=delazm+360.
          delazm=abs(delazm)
          IF(delazm < dazmin) THEN
            jj=j
            dazmin=delazm
          END IF
        END DO
        DO i=1,maxvelgate
          datarz(i,k) = max(velvol(i,jj,k),-30.)
        END DO
      END DO
 
      CALL xcolfil(datarz,rplt,zplt,iwrz,xvw,yvw, &
                   maxvelgate,maxvelgate,nvel_tilts,cl,ncl,mode)

      tem0 = azimuth_vel(iazimuth)+180.0
      IF( tem0 > 360.) tem0=tem0-360.

      DO k=1,nvel_tilts
        jj=1
        dazmin=720.
        DO j=1,maxazim
          delazm=azmrvol(j,k)-tem0
          IF(delazm > 180.) delazm=delazm-360.
          IF(delazm < -180.) delazm=delazm+360.
          delazm=abs(delazm)
          IF(delazm < dazmin) THEN
            jj=j
            dazmin=delazm
          END IF
        END DO
        DO i=1,maxvelgate
          datarz(i,k) = max(velvol(i,jj,k),-30.)
        END DO
      END DO
      CALL xcolfil(datarz,rpltm,zplt,iwrz,xvw,yvw,                   &
                   maxvelgate,maxvelgate,nvel_tilts,cl,ncl,mode)

      CALL xwdwof

      CALL xaxsca(-h_range_max, h_range_max, 0.1*h_range_max,         &
                       0.0, z_range_max, 0.5 )

      CALL xchsiz(0.05*z_range_max )  ! Set character size
      CALL xcharc((-h_range_max_w+h_range_max_e)*0.5, &
               h_range_max_n+2*h_range_max*0.07, trim(files(ifiles)) )

      WRITE(string,'(a,f6.2)')                                        &
            ' Radial Velocity at Azimuth Angle=',azimuth_vel(iazimuth)
      CALL xcharc(0.0, z_range_max*1.10, trim(string) )

      WRITE(stringtime,'(a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') 'First Scan Time:',&
            imon,'/',iday,'/',iyear,               &
            ' ',ihour,':',imin
      CALL xcharc(0.0, z_range_max*1.04, trim(stringtime) )

      CALL xchsiz(0.04*z_range_max )  ! Set character size
      CALL xcpalet(1)            ! Plot horizontal color palette

      CALL xframe

    END DO ! iazimuth

    DEALLOCATE( rplt )
    DEALLOCATE( rpltm )
    DEALLOCATE( zplt )
    DEALLOCATE( iwrz )
    DEALLOCATE( xvw )
    DEALLOCATE( yvw )
    DEALLOCATE( datarz )

    END IF   ! n_azimuth_vel
!
!-----------------------------------------------------------------------
! End of processing for this file
!-----------------------------------------------------------------------
!
    DEALLOCATE( elvref )
    DEALLOCATE( elvvel )
    DEALLOCATE( rngrvol )
    DEALLOCATE( azmrvol )
    DEALLOCATE( elvrvol )
    DEALLOCATE( refvol )
        !write(6,*) 'deallocating refvol'
    IF (ALLOCATED(vnyquist)) DEALLOCATE( vnyquist )
    DEALLOCATE( rngvvol )
    DEALLOCATE( azmvvol )
    DEALLOCATE( elvvvol )
    DEALLOCATE( velvol )

    IF(wrtuaref /= 0 .OR. wrtuavel /= 0) DEALLOCATE(uavol)
    IF(wrtvort  /= 0 ) DEALLOCATE( vortvol )
    IF(wrtdualp > 0) THEN
      DEALLOCATE( zdrvol )
      DEALLOCATE( kdpvol )
      DEALLOCATE( rhohvvol )
    END IF
    IF(wrtqx > 0) THEN
      DEALLOCATE( qrvol )
      DEALLOCATE( qsvol )
      DEALLOCATE( qhvol )
    END IF
    IF(wrtuvwt > 0) THEN
      DEALLOCATE( uvol )
      DEALLOCATE( vvol )
      DEALLOCATE( wvol )
      DEALLOCATE( tvol )
    END IF

  END DO ! ifiles
  IF( zxplot_initialized == 1 ) CALL xgrend  ! End graphics

  STOP  999

  100 CONTINUE
  Write(6,'(a)') 'Error reading namelist input file.'
  Write(6,'(a)') 'Program aborted.'
  STOP  101

END PROGRAM PLT_RADAR_TILT

!
!##################################################################
!##################################################################
!######                                                      ######
!######      Advanced Regional Prediction System (ARPS)      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE OUTPUT_REGULAR(regular_file, ntilts_max,ntilts,              &
                  radar_symbol,radar_lat,radar_lon,radar_elevation,     &
                  ngate,nray,fstgat,gatsp,int_razim,eleva,              &
                  volume )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Out put 3d regular polar grid data for further retrival.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Hu
!  03/19/02.
!
!  MODIFICATION HISTORY:
!
! INPUT:
!
!   regular_file  file name that data will be saved into
!   ntilts_max    maximum number of tilts
!   ntilts        actual number of tilts in this observation data
!   radar_symbol  radar Name 
!   radar_lat     radar location: latitude
!   radar_lon     radar location: longitude
!   radar_elevation  radar location: elavation
!
!   ngate      number of gates of observation data per radial
!   nray       number of radials of data( here = 361 )
!   fstgat     distance from radar of first gate of data(meters)
!   gatsp      gate spacing of data (meters)
!
!   volume(ngate,nray,ntilts)    Array containing one full volume scan of
!                                observation(reflectivity or radial velocity
!                                which has been intepolated to 1-degree
!                                regular spaced polar angles
!  int_razim(nray)     azimuth angle for each radial (degrees)
!  eleva(ntilts_max)   elevation angle for each tilt
!
!
!
! OUTPUT:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  INPUT RADAR PARAMETER AND OBSERVATION
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: regular_file
  INTEGER :: ntilts_max
  INTEGER :: ntilts
  CHARACTER(LEN=4) :: radar_symbol
  REAL :: radar_lat
  REAL :: radar_lon
  REAL :: radar_elevation

  INTEGER :: ngate
  INTEGER :: nray
  REAL ::    fstgat
  REAL ::    gatsp

  REAL :: volume(ngate,nray,ntilts)      ! Array containing one
                                         ! full volume scan of
                                         ! reflectivity intepolated to
                                         ! 1-degree spaced polar angles
  REAL :: int_razim(nray)
  REAL :: eleva(ntilts_max)

!
!-----------------------------------------------------------------------
!
!   WRITE(*,*) regular_file
!   WRITE(*,*) ntilts_max,ntilts,ngate,nray,ntilts
!   WRITE(*,*) radar_symbol,radar_lat,radar_lon,radar_elevation
!   WRITE(*,*) ngate,nray,fstgat,gatsp
!   WRITE(*,*) int_razim
!   WRITE(*,*) eleva
!  WRITE(*,*) volume
!
  OPEN(27,file=trim(regular_file),status='unknown',form='unformatted')
  WRITE(27) ntilts_max,ntilts
  WRITE(27) radar_symbol,radar_lat,radar_lon,radar_elevation
  WRITE(27) ngate,nray,fstgat,gatsp
  WRITE(27) int_razim
  WRITE(27) eleva
  WRITE(27) volume
  CLOSE(27)
  
  WRITE(*,*) ' The file:',trim(regular_file),                             &
                 'has been written out successfully'

END SUBROUTINE OUTPUT_REGULAR

SUBROUTINE pltradxs(ngate,nazim,ntilt,azimvol,varvol,            &
           rplt,rpltm,zplt,datarz,                               &
           varname,azimuth,rhimode,missing_value,                &
           xmin,xmax,xtick,zmin,zmax,ztick,                      &
           coltab,ibgncol,iendcol,coltabfn,pltmin,pltmax,        &
           iyear,imon,iday,ihour,imin,                           &
           cl,iwrz,xrw,yrw)
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: ngate,nazim,ntilt
  REAL, INTENT(IN) :: azimvol(nazim,ntilt)
  REAL, INTENT(IN) :: varvol(ngate,nazim,ntilt)
  REAL, INTENT(IN) :: rplt(ngate,ntilt)
  REAL, INTENT(IN) :: rpltm(ngate,ntilt)
  REAL, INTENT(IN) :: zplt(ngate,ntilt)
  REAL, INTENT(INOUT) :: datarz(ngate,ntilt)
  CHARACTER(LEN=30), INTENT(IN) :: varname
  REAL, INTENT(IN) :: azimuth
  INTEGER, INTENT(IN) :: rhimode
  REAL, INTENT(IN) :: missing_value
  REAL, INTENT(IN) :: xmin,xmax,xtick,zmin,zmax,ztick
  INTEGER, INTENT(IN) :: coltab
  INTEGER, INTENT(IN) :: ibgncol,iendcol
  CHARACTER(LEN=256), INTENT(IN) :: coltabfn
  REAL, INTENT(IN) :: pltmin,pltmax
  INTEGER, INTENT(IN) :: iyear,imon,iday,ihour,imin
  REAL, INTENT(INOUT) :: cl(100)
  INTEGER, INTENT(INOUT) :: iwrz(ngate,ntilt)
  REAL, INTENT(INOUT) :: xrw(8*ngate)
  REAL, INTENT(INOUT) :: yrw(8*ngate)
!
! Misc local variables
!
  INTEGER :: ncl,mode
  INTEGER :: i,j,jj,k
  REAL :: dazmin,delazm
  REAL :: tem0,xlab
  REAL :: varmin,varmax
  REAL :: pltminlim,pltmaxlim
  CHARACTER(LEN=128) :: string 
  CHARACTER(LEN=128) :: stringtime

  varmin=999.
  varmax=-999.
  pltminlim=pltmin-(cl(2)-cl(1))
  pltmaxlim=pltmax+(cl(2)-cl(1))

  CALL xpspac(0.1,0.9,0.2,0.6)
  CALL xmap(xmin,xmax,zmin,zmax)
  CALL xlabmask(1)

  print *, 'XS RHI mode: ',rhimode
  print *, 'XS Color table is: ',TRIM(coltabfn)
  IF( coltab == -1 ) CALL xstctfn(coltabfn)
  CALL xsetclrs(coltab)
  CALL xcolor(1)    ! Set color to black
  CALL xcfont(3)
!
! Set tentative color interval and limit on the min/max number of contours
!
  CALL xnctrs(3,50)  ! Set lower and upper limits
                            ! on the number of contours allowed
                            ! It's used by both xcolfil and xclimit
  ncl=1
  mode=2
  xlab=0.5*(xmin+xmax)

  CALL xctrbadv(1)   ! Turn on bad value checking
                     ! in contouring routine
  CALL xbadval( missing_value ) ! Specify bad value flag
!
! plot a color filled contour field
!
  print *, ' XS color indices: ',ibgncol,iendcol
  CALL xctrclr(ibgncol,iendcol)

  print *, ' XS ctr limits: ',pltmin,pltmax
  CALL xctrlim(pltmin,pltmax) 
  CALL xwindw(xmin,xmax,zmin,zmax)

  datarz=missing_value

  print *, '  XS azimuth: ',azimuth
  DO k=1,ntilt
    jj=1
    dazmin=720.
    DO j=1,nazim
      delazm=azimvol(j,k)-azimuth
      IF(delazm > 180.) delazm=delazm-360.
      IF(delazm < -180.) delazm=delazm+360.
      delazm=abs(delazm)
      IF(delazm < dazmin) THEN
        jj=j
        dazmin=delazm
      END IF
    END DO
    print *, ' elev k=',k,'  nearest azimuth: ',azimvol(jj,k)
    DO i=1,ngate
      IF(varvol(i,jj,k) > -999.) THEN
        varmin=min(varmin,varvol(i,jj,k))
        varmax=max(varmax,varvol(i,jj,k))
        datarz(i,k) = min(max(varvol(i,jj,k),pltminlim),pltmaxlim)
      END IF
    END DO
  END DO
  
  print *, ' Here XS 0'
  print *, ' ngate: ',ngate
  print *, ' ntilt: ',ntilt
  print *, ' ncl: ',ncl
  print *, ' mode: ',mode
 
  CALL xcolfil(datarz,rplt,zplt,iwrz,xrw,yrw, &
                   ngate,ngate,ntilt,cl,ncl,mode)
  print *, ' Here XS 1'

  IF(rhimode < 1) THEN
    tem0 = azimuth+180.0
    IF( tem0 > 360.) tem0=tem0-360.
    print *, '  XS azimuth 2: ',tem0

    datarz=missing_value

    DO k=1,ntilt
      jj=1
      dazmin=720.
      DO j=1,nazim
        delazm=azimvol(j,k)-tem0
        IF(delazm > 180.) delazm=delazm-360.
        IF(delazm < -180.) delazm=delazm+360.
        delazm=abs(delazm)
        IF(delazm < dazmin) THEN
          jj=j
          dazmin=delazm
        END IF
      END DO
      print *, ' elev k=',k,'  nearest azimuth: ',azimvol(jj,k)
      DO i=1,ngate
        IF(varvol(i,jj,k) > -999.) THEN
          varmin=min(varmin,varvol(i,jj,k))
          varmax=max(varmax,varvol(i,jj,k))
          datarz(i,k) = min(max(varvol(i,jj,k),pltminlim),pltmaxlim)
        END IF
      END DO
    END DO

    print *, ' Here XS 3'
    CALL xcolfil(datarz,rpltm,zplt,iwrz,xrw,yrw,           &
                   ngate,ngate,ntilt,cl,ncl,mode)
    print *, ' Here XS 4'
  END IF

  CALL xwdwof

  CALL xaxsca(xmin,xmax,xtick,zmin,zmax,ztick)

  WRITE(string,'(a,a,f6.2)') TRIM(varname),                &
            ' at azimuth angle=',azimuth
  CALL xchsiz(0.06*zmax )  ! Set character size
  CALL xcharc(xlab,zmax*1.10, trim(string) )

  WRITE(stringtime, &
   '(i2.2,a,i2.2,a,i4,1x,i2.2,a,i2.2,4x,a,f7.1,a,f7.1)')   &
     imon,'/',iday,'/',iyear,ihour,':',imin, &
     '  min:',varmin,'  max:',varmax
  WRITE(6,'(a,f7.1,a,f7.1)')   '  min:',varmin,'  max:',varmax
  CALL xchsiz(0.05*zmax )  ! Set character size
  CALL xcharc(xlab, zmax*1.04, trim(stringtime) )

  CALL xchsiz(0.04*zmax )  ! Set character size
  CALL xcpalet(1)            ! Plot horizontal color palette

  CALL xframe
  print *, ' Here XS 5'
  RETURN
END SUBROUTINE pltradxs
!
SUBROUTINE beamhgt(elvang,range,height,sfcrng)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the height of the radar beam and the along-
!  ground distance from the radar as a function
!  distance along the radar beam (range) and radar
!  elevation angle (elvang).
!
!  This method assumes dn/dh is constant such that the
!  beam curves with a radius of 4/3 of the earth's radius.
!  This is from Eq. 2.28 of Doviak and Zrnic', Doppler Radar
!  and Weather Observations, 1st Ed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  06/22/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    elvang   Elevation angle (degrees) of radar beam
!    range    Distance (meters) along radar beam from radar
!
!  OUTPUT:
!    height   Height (meters) of beam above ground.
!    sfcrng   Distance (meters) of point along ground from radar.
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
  REAL, INTENT(IN) :: elvang
  REAL, INTENT(IN) :: range
  REAL, INTENT(OUT) :: height
  REAL, INTENT(OUT) :: sfcrng
!
  DOUBLE PRECISION :: eradius,frthrde,eighthre,fthsq,deg2rad
  PARAMETER (eradius=6371000.,                                          &
             frthrde=(4.*eradius/3.),                                   &
             eighthre=(8.*eradius/3.),                                  &
             fthsq=(frthrde*frthrde),                                   &
             deg2rad=(3.14592654/180.))
!
  DOUBLE PRECISION :: elvrad,hgtdb,rngdb,drange
!
  elvrad=deg2rad*DBLE(elvang)
  drange=DBLE(range)
  hgtdb = SQRT(drange*drange + fthsq +                                  &
                eighthre*drange*SIN(elvrad)) -                          &
                frthrde
  height=hgtdb
  rngdb = frthrde *                                                     &
           ASIN (drange*COS(elvrad)/(frthrde + hgtdb) )
  sfcrng=rngdb
  RETURN
END SUBROUTINE beamhgt
