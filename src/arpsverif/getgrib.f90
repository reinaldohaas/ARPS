!Edited to use FORTRAN 77 include files...EMK
!arpscntl RLC istatus
!########################################################################
!########################################################################
!######                                                            ######
!######                    SUBROUTINE GETGRIB                      ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE GETGRIB(nx_ext,ny_ext,nz_ext,                                &
                        dir_extd,extdname,extdopt,                      &
                        extdfmt,                                        &
                        extdinit,extdfcst,julfname,                     &
                        iproj_ext,scale_ext,                            &
                        trlon_ext,latnot_ext,x0_ext,y0_ext,             &
                        lat_ext,lon_ext,                                &
                        p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,         &
                        qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,             &
                        tsfc_ext,tsoil_ext,                             &
                        wetsfc_ext,wetdp_ext,wetcanp_ext,               &
                        snowcvr_ext,trn_ext,psfc_ext,                   &
                        T_2m_ext,RH_2m_ext,U_10m_ext,                   &
                        V_10m_ext,MSLP_ext,RH_ext,                      &
                        undrgrnd,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads in a NMC GRIB file (RUC #87, RUC #211, Eta #212, or RUC2 #236)
!  for processing by CVT2VERIF.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Eric Kemp
!  November 1999
!  Based on subroutines GETNMCRUC87, GETNMCRUC211, GETNMCETA212,
!  GETNMCRUCN236, and GETNMCRUCP236.
!
!  MODIFICATION HISTORY:
!  Eric Kemp, 14 December 1999
!  Added correction for converting relative humidity to specific
!  humidity.  Relative humidity is now always treated with respect
!  to liquid water.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext,ny_ext,nz_ext  The dimension of data arrays
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format

!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsfc_ext      Surface temperature
!    tsoil_ext     Soil temperature
!    wetsfc_ext    Top layer soil moisture
!    wetdp_ext     Deep soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    istatus       status indicator (0=success, -1=unknown error,
!		   1=file not found, 2=data error)
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!
!-----------------------------------------------------------------------
!
!  Use modules.
!
!-----------------------------------------------------------------------

!  USE gribcst2
!  USE globcst
!  USE phycst

  USE gribcst2

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'

  INTEGER :: nx_ext,ny_ext,nz_ext

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname

!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------

  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext

!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL :: tsfc_ext   (nx_ext,ny_ext)      ! Temperature at surface (K)
  REAL :: tsoil_ext  (nx_ext,ny_ext)      ! Deep soil temperature (K)
  REAL :: wetsfc_ext (nx_ext,ny_ext)      ! Surface soil moisture
  REAL :: wetdp_ext  (nx_ext,ny_ext)      ! Deep soil moisture
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  INTEGER :: snowcvr_ext(nx_ext,ny_ext)   ! Snow cover

  REAL :: trn_ext    (nx_ext,ny_ext)      ! External terrain (m)
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

  REAL :: T_2m_ext(nx_ext,ny_ext)
  REAL :: RH_2m_ext(nx_ext,ny_ext)
  REAL :: U_10m_ext(nx_ext,ny_ext)
  REAL :: V_10m_ext(nx_ext,ny_ext)
  REAL :: MSLP_ext(nx_ext,ny_ext)
  REAL :: RH_ext(nx_ext,ny_ext,nz_ext)

!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------

  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)

  INTEGER :: istatus

!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------

  INTEGER :: iproj
  REAL :: scale,trlon,x0,y0
  REAL :: latnot(2)

!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=80) :: gribfile
  CHARACTER (LEN=13) :: gribtime
  INTEGER :: i,j,k,kk
  INTEGER :: iyr,imo,iday,myr, jldy
  INTEGER :: ihr,imin,isec
  INTEGER :: ifhr,ifmin,ifsec
  INTEGER :: grbflen, len1, lenrun

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  REAL :: govrd

  INTEGER :: chklev, lvscan

  INTEGER :: iret             ! Return flag

  LOGICAL :: exist

!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
!  REAL :: zk_grb(nzgrb)   ! Vertical coordinate parameters
  REAL :: zk_grb(nz_ext)   ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
!  INTEGER :: nit(nzgrb)   ! Number of x-points for thinned grid
  INTEGER :: nit(nz_ext)   ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth

  INTEGER :: var_lev3d(nz_ext,n3dvs,n3dlvt)
  REAL :: rcdata(nx_ext*ny_ext)      ! temporary data array
  REAL :: var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt) ! GRIB variables
  REAL :: var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt) ! GRIB 3-D
                                                       ! variables
  LOGICAL :: ibms(nx_ext*ny_ext)  ! BMS logical array for data bit map

  REAL :: landseamask(nx_ext,ny_ext)
  REAL :: theta_ext(nx_ext,ny_ext,nz_ext)
  REAL :: msf_ext(nx_ext,ny_ext,nz_ext) ! Montgomery stream function
  REAL :: vpt_ext(nx_ext,ny_ext,nz_ext) ! Virtual potential temperature
  REAL :: mixrat_ext(nx_ext,ny_ext,nz_ext) ! Water vapor mixing ratio
  REAL :: cplp(nx_ext,ny_ext,nz_ext) ! condensation pressure of lifted
                                     ! parcel (Pa)
  REAL :: tv_ext, tvc_ext
  REAL :: ROVCP_P, CPD_P, G0_P, RD_P
  REAL :: tema, temb
  INTEGER :: ii,jj
  REAL :: A,B
  REAL :: qvsat, pilev, qvsatice

  INTEGER :: undrgrnd(nx_ext,ny_ext,nz_ext)

!-----------------------------------------------------------------------
!
! Functions
!
!-----------------------------------------------------------------------

  REAL,EXTERNAL :: f_qvsatl
!fpp$ expand (f_qvsatl)
!dir$ inline always f_qvsatl

!-----------------------------------------------------------------------
!
! Parameters
!
!-----------------------------------------------------------------------

  REAL,PARAMETER :: missing = -9999.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = -1

  DO j = 1,n3dlvt
    DO i = 1,n3dvs
      var_id3d(i,j) = 255 ! Initialized as undefined.  EMK
    END DO
  END DO

  DO j = 1,n2dlvt
    DO i = 1,n2dvs
      var_id2d(i,j) = 255 ! Initialized as undefined.  EMK
    END DO
  END DO

  IF(extdfcst == '         ') extdfcst='000:00:00'

  READ (extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                  &
       iyr,imo,iday,ihr,imin,isec

  CALL julday(iyr,imo,iday,jldy)

  myr=MOD(iyr,100)
  ifhr=0
  ifmin=0
  ifsec=0

  READ(extdfcst,'(i3)',ERR=4,END=4) ifhr

  4 CONTINUE

  WRITE (gribtime,'(i4.4,i2.2,i2.2,i2.2,a1,i2.2)') iyr,imo,iday,ihr,'f',ifhr

  len1=LEN(dir_extd)
  grbflen=len1

  CALL strlnth( dir_extd, grbflen )

  IF( grbflen == 0 .OR. dir_extd(1:grbflen) == ' ' ) THEN
    dir_extd = '.'
    grbflen=1
  END IF

  IF( dir_extd(grbflen:grbflen) /= '/'.AND.grbflen < len1 ) THEN
    grbflen=grbflen+1
    dir_extd(grbflen:grbflen)='/'
  END IF

  lenrun = LEN( extdname )
  CALL strlnth( extdname, lenrun )

  gribfile = dir_extd(1:grbflen)//extdname(1:lenrun)                    &
                                //'.'//gribtime(1:13)
  grbflen = grbflen + lenrun + 14

  INQUIRE (FILE=gribfile, EXIST=exist)
  IF (.NOT. exist) THEN
    PRINT *, 'getgrib: ERROR: File does not exist: ', TRIM(gribfile)
    istatus = 1
    RETURN
  END IF

!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------

  govrd = g/rd

! EMK:  This could perhaps be put in a separate subroutine...

  IF (extdopt == 1) THEN
    gridtyp  = ruc87grid
    mproj_grb = ruc87proj

    n2dvars  = ruc87nvs2d
    n2dlvtps = ruc87nlvt2d

    DO k=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,k) = ruc87var_id2d(n,k)
      END DO
      levtyp2d(k) = ruc87levs2d(k)
    END DO

    n3dvars  = ruc87nvs3d
    n3dlvtps = ruc87nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = ruc87var_id3d(n,m)
      END DO
      levtyp3d(m) = ruc87levs3d(m)
    END DO

  ELSE IF (extdopt == 2) THEN
    gridtyp  = eta212grid
    mproj_grb = eta212proj

    n2dvars  = eta212nvs2d
    n2dlvtps = eta212nlvt2d

    DO k=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,k) = eta212var_id2d(n,k)
      END DO
      levtyp2d(k) = eta212levs2d(k)
    END DO

    n3dvars  = eta212nvs3d
    n3dlvtps = eta212nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = eta212var_id3d(n,m)
      END DO
      levtyp3d(m) = eta212levs3d(m)
    END DO
  ELSE IF (extdopt == 7) THEN
    gridtyp  = ruc211grid
    mproj_grb = ruc211proj

    n2dvars  = ruc211nvs2d
    n2dlvtps = ruc211nlvt2d

    DO k=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,k) = ruc211var_id2d(n,k)
      END DO
      levtyp2d(k) = ruc211levs2d(k)
    END DO

    n3dvars  = ruc211nvs3d
    n3dlvtps = ruc211nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = ruc211var_id3d(n,m)
      END DO
      levtyp3d(m) = ruc211levs3d(m)
    END DO

  ELSE IF (extdopt == 11) THEN

    gridtyp  = rucn236grid
    mproj_grb = rucn236proj

    n2dvars  = rucn236nvs2d
    n2dlvtps = rucn236nlvt2d

    DO k=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,k) = rucn236var_id2d(n,k)
      END DO
      levtyp2d(k) = rucn236levs2d(k)
    END DO

    n3dvars  = rucn236nvs3d
    n3dlvtps = rucn236nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = rucn236var_id3d(n,m)
      END DO
      levtyp3d(m) = rucn236levs3d(m)
    END DO

  ELSE IF (extdopt == 12) THEN

    gridtyp  = rucp236grid
    mproj_grb = rucp236proj

    n2dvars  = rucp236nvs2d
    n2dlvtps = rucp236nlvt2d

    DO k=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,k) = rucp236var_id2d(n,k)
      END DO
      levtyp2d(k) = rucp236levs2d(k)
    END DO

    n3dvars  = rucp236nvs3d
    n3dlvtps = rucp236nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = rucp236var_id3d(n,m)
      END DO
      levtyp3d(m) = rucp236levs3d(m)
    END DO

  ENDIF

!  CALL rdnmcgrb(gribfile,grbflen, gribtime,                             &
!                gridesc, iproj_grb, gthin,                              &
!                ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
!                pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
!                latsw,lonsw, latne,lonne,                               &
!                latrot,lonrot,angrot,                                   &
!                latstr,lonstr,facstr,                                   &
!                lattru1,lattru2,lontrue,                                &
!                scanmode, iscan,jscan,kscan,                            &
!                ires,iearth,icomp,                                      &
!                jpenta,kpenta,mpenta,ispect,icoeff,                     &
!                xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
!                iret)

  CALL rdnmcgrb2(gribfile,grbflen, gribtime,                            &
                gridesc, iproj_grb, gthin,                              &
                ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                latsw,lonsw, latne,lonne,                               &
                latrot,lonrot,angrot,                                   &
                latstr,lonstr,facstr,                                   &
                lattru1,lattru2,lontrue,                                &
                scanmode, iscan,jscan,kscan,                            &
                ires,iearth,icomp,                                      &
                jpenta,kpenta,mpenta,ispect,icoeff,                     &
                xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                nx_ext,ny_ext,nz_ext,var_lev3d,var_grb3d,               &
                var_grb2d,ibms,rcdata,iret)

  max_nr2d = 0
  DO n=1,n2dvars
    DO m=1,n2dlvtps
      max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
    END DO
  END DO

  max_nr3d = 0
  DO n=1,n3dvars
    max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
  END DO

  IF ( max_nr3d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 3-D variable was found in the GRIB file',                   &
        'Program stopped in GETGRIB.'
    istatus = 2
    RETURN
  END IF

  IF ( max_nr2d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 2-D variables was found in the GRIB file'
  END IF

!arpscntl RLC print
  WRITE (6,'(/a7,2x,6(i7))')                                            &
      'Lev\VID',(var_id3d(n,1),n=1,n3dvars)
  DO k=1,max_nr3d
    WRITE (6,'(i5,4x,6(i7))')                                          &
        k,(var_lev3d(k,n,1),n=1,n3dvars)
  END DO

  DO k=1,max_nr3d
    DO n=2,n3dvars
      IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
        WRITE (6,'(a)')                                                 &
            'Variables were not at the same level.',                    &
            'Program stopped in GETGRIB.'
        istatus = 2
	RETURN
      END IF
    END DO
  END DO

  IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
    iproj_ext = 1
  ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -1
  ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
    iproj_ext = 2
  ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -2
  ELSE IF ( iproj_grb == 1 ) THEN
    iproj_ext = 3
  ELSE IF ( iproj_grb == 0 ) THEN
    iproj_ext = 4
  ELSE
    WRITE (6,'(a)')                                                     &
        'Unknown map projection. Set to non-projection.'
    iproj_ext = 0
  END IF

  scale_ext = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)

!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------

! EMK:  This could perhaps be put into a separate subroutine...

  WRITE(6,*)'Retrieving 2-D variables...'
  DO j=1,ny_ext
    DO i=1,nx_ext
      landseamask(i,j) = real(0)
      psfc_ext(i,j) = missing
      trn_ext(i,j) = missing
      tsfc_ext(i,j) = missing
      T_2m_ext(i,j) = missing
      RH_2m_ext(i,j) = missing
      U_10m_ext(i,j) = missing
      V_10m_ext(i,j) = missing
      snowcvr_ext(i,j) = -9999
      wetcanp_ext(i,j) = missing
      tsfc_ext(i,j) = missing
      tsoil_ext(i,j) = missing
      wetsfc_ext(i,j) = missing
      wetdp_ext(i,j) = missing
    END DO
  END DO

  DO j=1,ny_ext
    DO i=1,nx_ext
      DO jj = 1,n2dlvt
      DO ii = 1,n2dvs
        IF (var_id2d(ii,jj) == 81) THEN ! Land/sea mask
          landseamask(i,j) = var_grb2d(i,j,ii,jj)
        ENDIF
      END DO
      END DO
    END DO
  END DO

  DO j=1,ny_ext
    DO i=1,nx_ext

      DO jj = 1,n2dlvt
      DO ii = 1,n2dvs

      SELECT CASE (var_id2d(ii,jj))
      CASE (1) ! Pressure (Pa)
        IF ( levtyp2d(jj) == 001 ) THEN ! Surface pressure (Pa)
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            psfc_ext   (i,j) = var_grb2d(i,j,ii,jj)
          END IF
        ENDIF
      CASE (2) ! Pressure reduced to MSL (Pa)
        IF ( levtyp2d(jj) == 102) THEN
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            MSLP_ext(i,j) = var_grb2d(i,j,ii,jj)/100.0 ! Pa to mb
          ENDIF
        ENDIF
      CASE (7) ! Geopotential height (gpm)
        IF (levtyp2d(jj) == 001 ) THEN ! Surface geopotential height
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            trn_ext    (i,j) = var_grb2d(i,j,ii,jj)
          ENDIF
        ENDIF
      CASE (8) ! Geometric height (m)
        IF (levtyp2d(jj) == 001 ) THEN ! Surface geometric height
          IF ( var_nr2d(ii,jj) /= 0) THEN
            trn_ext    (i,j) = var_grb2d(i,j,ii,jj)
          ENDIF
        ENDIF
      CASE (11) ! Temperature (K)
        IF (levtyp2d(jj) == 105 ) THEN ! 2 m temperature (K)
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            T_2m_ext(i,j)= var_grb2d(i,j,ii,jj)
          ENDIF
        ELSE IF (levtyp2d(jj) == 001 ) THEN ! Ground temperature (K)
          tsfc_ext(i,j) = var_grb2d(i,j,ii,jj)
          tsoil_ext(i,j) = var_grb2d(i,j,ii,jj)
        ENDIF
      CASE (33) ! U-winds (m/s)
        IF (levtyp2d(jj) == 105 ) THEN ! 10 m U-winds (m/s)
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            U_10m_ext(i,j)= var_grb2d(i,j,ii,jj)
          ENDIF
        ENDIF
      CASE (34) ! V-winds (m/s)
        IF (levtyp2d(jj) == 105 ) THEN ! 10 m V-winds (m/s)
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            V_10m_ext(i,j)= var_grb2d(i,j,ii,jj)
          ENDIF
        ENDIF
      CASE (52) ! Relative Humidity (%)
        IF (levtyp2d(jj) == 105 ) THEN ! 2 m relative humidity (%)
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            rh_2m_ext(i,j)= var_grb2d(i,j,ii,jj)
          ENDIF
        ENDIF
      CASE (66) ! Snow depth (m)
        IF (levtyp2d(jj) == 001 ) THEN
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            IF (var_grb2d(i,j,ii,jj) == 2.54e-4) then ! >= 0.01 inches
              snowcvr_ext(i,j) = 1
            ELSE
              snowcvr_ext(i,j) = 0
            ENDIF
          ENDIF
        ENDIF
      CASE (129) ! MAPS mean sea level pressure (Pa)
        IF (levtyp2d(jj) == 102) THEN
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            MSLP_ext(i,j) = var_grb2d(i,j,ii,jj)/100.0 ! Pa to mb
          ENDIF
        ENDIF
      CASE (130) ! Eta mean sea level pressure (Pa)
        IF (levtyp2d(jj) == 102) THEN
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            MSLP_ext(i,j) = var_grb2d(i,j,ii,jj)/100.0 ! Pa to mb
          ENDIF
        ENDIF
      CASE (143) ! Categorical snow (yes=1;no=0)
        IF (levtyp2d(jj) == 001 ) THEN
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            snowcvr_ext(i,j) = nint(var_grb2d(i,j,ii,jj))
          ENDIF
        ENDIF
      CASE (223) !  Plant canopy surface water (kg/m**2)
        IF (levtyp2d(jj) == 001 ) THEN
          IF ( var_nr2d(ii,jj) /= 0 ) THEN
            wetcanp_ext(i,j) = var_grb2d(i,j,4,1)*1.e-3     ! in meter
          ENDIF
        ENDIF

      CASE DEFAULT
      END SELECT

      END DO
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------

  CPD_P   = 1004.686    ! cp in RUC
  ROVCP_P = 0.285714    ! rd/cp used in RUC
  G0_P    = 9.80665     ! gravity in RUC

  nz1 = MIN(var_nr3d(1,1),nz_ext)

  IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
    chklev = 1
    lvscan = 0
  ELSE
    chklev = -1
    lvscan = nz1+1
  END IF

  WRITE(6,*)'Retrieving 3-D variables...'
  DO k = 1,nz_ext
    DO j = 1,ny_ext
      DO i = 1,nx_ext
        p_ext(i,j,k) = missing
        hgt_ext(i,j,k) = missing
        u_ext(i,j,k) = missing
        v_ext(i,j,k) = missing
        qv_ext(i,j,k) = missing
        qc_ext(i,j,k) = missing
        qr_ext(i,j,k) = missing
        qi_ext(i,j,k) = missing
        qs_ext(i,j,k) = missing
        qh_ext(i,j,k) = missing
        msf_ext(i,j,k) = missing
        theta_ext(i,j,k) = missing
        vpt_ext(i,j,k) = missing
        undrgrnd(i,j,k) = -999
      END DO
    END DO
  END DO

  DO k = 1,nz1
    kk = chklev * k + lvscan
    DO j=1,ny_ext
      DO i=1,nx_ext
        DO jj = 1,n3dlvt
        DO ii = 1,n3dvs
          SELECT CASE (var_id3d(ii,jj))
          CASE (1) ! Pressure (Pa)
            IF (levtyp3d(jj) == 109) THEN ! Pressure at Hybrid Sfc.
              p_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
            ENDIF
          CASE (7) ! Height (m)
            hgt_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)    ! height (m)
          CASE (11) ! Temperature
            t_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (13) ! Potential temperature (K)
            theta_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (33) ! u wind (m/s)
            u_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (34) ! v wind (m/s)
            v_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (37) ! Montgomery stream function
            msf_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (51) ! Specific humidity
            qv_ext (i,j,kk) = var_grb3d(i,j,k,ii,jj)    ! Spec. humidity
                                                        ! (kg/kg)
          CASE (52) ! Relative humidity (%)
            IF (levtyp3d(jj) == 100) THEN ! pressure surface
              RH_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
            ENDIF
          CASE (85) ! Soil temperature (K)
            IF (levtyp3d(jj) == 111) THEN ! soil data
              IF (extdopt == 1) THEN ! RUC #87
                IF (var_nr3d(ii,jj) /= 0) THEN
                  tsfc_ext   (i,j) = var_grb3d(i,j,1,ii,jj)
                  tsoil_ext  (i,j) = var_grb3d(i,j,3,ii,jj)
                ENDIF
              ELSE IF (extdopt == 2) THEN ! Eta #212
                IF (var_nr3d(ii,jj) /= 0) THEN
                  IF ( nint(landseamask(i,j)).eq.1 ) THEN
                           ! soil temp over land
                    tsoil_ext (i,j) = 0.1 * var_grb3d(i,j,1,ii,jj)       &
                            + 0.3 * var_grb3d(i,j,2,ii,jj)               &
                            + 0.6 * var_grb3d(i,j,3,ii,jj)
                            ! 0-100cm average
                  ENDIF
                ENDIF
              ELSE IF (extdopt == 11) THEN ! RUC Native #236
                IF (var_nr3d(ii,jj) /= 0) THEN
                  tsfc_ext   (i,j) = var_grb3d(i,j,1,ii,jj)
                  tsoil_ext  (i,j) = var_grb3d(i,j,3,ii,jj)
                ENDIF
              ENDIF
            ENDIF
          CASE (144) ! Volumetric soil moisture (m**3/m**3)
            IF (levtyp3d(jj) == 112) THEN ! soil data
              IF (extdopt == 1) THEN ! RUC # 87
                IF (var_nr3d(ii,jj) /= 0) THEN
                  wetsfc_ext (i,j) = var_grb3d(i,j,1,ii,jj)
                  wetdp_ext  (i,j) = var_grb3d(i,j,3,ii,jj)
                ENDIF
              ELSE IF (extdopt == 2) THEN ! Eta #212
                IF (var_nr3d(ii,jj) /= 0) THEN
                  wetsfc_ext(i,j) = var_grb3d(i,j,1,ii,jj)
                  wetdp_ext(i,j)  = 0.1 * var_grb3d(i,j,1,ii,jj)       &
                          + 0.3 * var_grb3d(i,j,2,ii,jj)               &
                          + 0.6 * var_grb3d(i,j,3,ii,jj) ! 0-100cm average
                ENDIF
              ELSE IF (extdopt == 11) THEN
                IF (var_nr3d(ii,jj) /= 0) THEN ! RUC Native #236
                  wetsfc_ext (i,j) = var_grb3d(i,j,1,ii,jj)
                  wetdp_ext  (i,j) = var_grb3d(i,j,3,ii,jj)
                ENDIF
              ENDIF
            ENDIF
          CASE (153) ! cloud water mixing ratio (kg/kg)
            qc_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (159) ! condensation pressure of lifted parcel (Pa)
            cplp(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (170) ! rain water mixing ratio (kg/kg)
            qr_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (171) ! snow mixing ratio (kg/kg)
            qs_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (178) ! ice mixing ratio (kg/kg)
            qi_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (179) ! graupel mixing ratio (kg/kg)
            qh_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (185) ! water vapor mixing ratio (kg/kg)
            mixrat_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE (189) ! virtual potential temperature (K)
            vpt_ext(i,j,kk) = var_grb3d(i,j,k,ii,jj)
          CASE DEFAULT
          END SELECT
        END DO
        END DO
      END DO
    END DO
  END DO

  WRITE(6,*)'Almost done with 3-D variables...'

!-----------------------------------------------------------------------
!
! If processing isobaric data, determine pressure at each constant
! pressure level and check for portions of constant pressure grids
! that are below the surface.
!
!-----------------------------------------------------------------------

  DO k=1,nz1
    kk = chklev * k + lvscan
    DO j=1,ny_ext
      DO i=1,nx_ext
        IF (extdopt == 2 .OR. extdopt == 12) THEN ! Pressure on constant
                                                  ! pressure surface.
          p_ext  (i,j,kk) = 100.0 * FLOAT(var_lev3d(k,1,1))
          IF (((p_ext(i,j,kk) > psfc_ext(i,j))) .AND.   &
              (psfc_ext(i,j) > 0.) ) THEN
            undrgrnd(i,j,kk) = 1
          ELSE
            undrgrnd(i,j,kk) = 0
          ENDIF
          IF( (extdopt == 12 .AND. p_ext(i,j,kk) > 0.0) .AND.        &
              (t_ext(i,j,kk) > 0.0) ) then
            qvsat= f_qvsatl( p_ext(i,j,kk), t_ext(i,j,kk) )
            qv_ext(i,j,kk)= RH_ext(i,j,kk) * qvsat * 0.01
          ENDIF
          IF( (extdopt == 2 .AND. p_ext(i,j,kk) > 0.0) .AND.        &
              (t_ext(i,j,kk) > 0.0) ) then
            qvsat= f_qvsatl( p_ext(i,j,kk), t_ext(i,j,kk) )
            RH_ext(i,j,kk) = qv_ext(i,j,kk)/(qvsat*0.01)
          ENDIF
        ENDIF
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
! If processing RUC native coordinate data, calculate temperature,
! height, and specific humidity.
!
!-----------------------------------------------------------------------

  DO k=1,nz1
    kk = chklev * k + lvscan
    DO j=1,ny_ext
      DO i=1,nx_ext
        IF (extdopt == 1) THEN
          pilev = (p_ext(i,j,kk)/100000.)**ROVCP_P
          tv_ext = theta_ext(i,j,kk)*pilev     ! Virtual Temperature (K)
          hgt_ext(i,j,kk) = (msf_ext(i,j,kk)-CPD_P*tv_ext)/G0_P
                                ! Height (m) from Montgomery function with
                                ! M = CpTv + gz

          tvc_ext = theta_ext(i,j,kk)                                   &
                       * (cplp(i,j,kk)/100000.)**ROVCP_P
                              ! Virtual Temperature (K) at LCL
          CALL tv2tq( tvc_ext, cplp(i,j,kk), tema, temb )
          t_ext(i,j,kk)  = tema                                         &
             * (p_ext(i,j,kk)/cplp(i,j,kk))**ROVCP_P ! Temperature (K)
          qv_ext(i,j,kk) = temb                      ! Specific humidity

        ENDIF
        IF (extdopt == 11) THEN
          A = real(100000)/p_ext(i,j,kk)
          A = A**ROVCP_P
          tvc_ext = vpt_ext(i,j,kk)/A ! Virtual Temperature
          B = 0.61*mixrat_ext(i,j,kk)
          B = real(1) + B
          t_ext(i,j,kk) = tvc_ext/B ! Temperature (K)

          A = mixrat_ext(i,j,kk)*p_ext(i,j,kk)
          A = A/(0.622 - mixrat_ext(i,j,kk))
          qv_ext(i,j,kk) = A*0.622/p_ext(i,j,kk) !SpecificHumidity
        ENDIF
      END DO
    END DO
  END DO

  WRITE(6,*)'Finished with 3-D variables...'

!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north, if necessary.
!
!-----------------------------------------------------------------------

  IF (extdopt == 1 .OR. extdopt == 2 .OR. extdopt == 7 .OR.             &
      extdopt == 11 .OR. extdopt == 12) THEN
    DO k=1,nz1
      CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),             &
                 lon_ext,u_ext(1,1,k),v_ext(1,1,k))
    END DO

!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------

    CALL setmapr(iproj,scale,latnot,trlon)
    CALL setorig(1,x0,y0)

  ENDIF

  istatus = 0

  RETURN
END SUBROUTINE getgrib

