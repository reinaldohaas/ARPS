!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDRADCOL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdradcol(nx,ny,nz,nsrc_rad,nvar_radin,                       &
           mx_rad,nz_rdr,mx_colrad,mx_pass,                             &
           raduvobs,radrhobs,radistride,radkstride,cloudopt,            &
           xs,ys,zs,latgr,longr,                                        &
           iuserad,iusechk,npass,nradfil,radfname,                      &
           srcrad,isrcrad,stnrad,latrad,lonrad,elvrad,                  &
           latradc,lonradc,irad,nlevrad,hgtradc,obsrad,                 &
           iprocv,jprocv,iradvel,indexrad,oindexrad,                    &
           refmax_mos_3d,ref_mos_3d,                                    &
           rhv_mos_3d,zdr_mos_3d,kdp_mos_3d,dist_mos_3d,                &
           ncolrad,ncolrad_mpi,istatus,istat_radar,                     &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads radar data stored as columns on a remapped grid.
!  Data stored in this way can be easily distributed to appropriate
!  processors in MPI algoritms.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  August, 1995
!
!  MODIFICATION HISTORY:
!
!  04/25/02 (Keith Brewster)
!  Added stride processing for data thinning, when desired.
!
!  05/04/02 (Leilei Wang and Keith Brewster)
!  Added reading of hdf formatted files.
!
!  10/29/2002 (Keith Brewster)
!  Improvement to alloc and error handling of hdf i/o.
!
!  2007 (Yunheng Wang)
!  Add support for new format radar file.
!
!  Nov 2007 (Kevin W. Thomas)
!  Suppress repeated messages in MPI mode.  Fix seg fault in MPI mode.
!
!  April, 2010 (Keith Brewster)
!  Added processing for reflectivity mosaic to avoid a re-reading files
!  later for the cloud analysis.
!
!  September, 2012 (Keith Brewster)
!  Added processing for reading and mosaiking of dual-pol variables.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nsrc_rad   number of sources of radar data
!    nvar_radin number of variables in the obsrad array
!    mx_rad     maximum number of radars
!    nz_rdr     maximum number of levels in a radar column
!    mx_colrad  maximum number of radar columns
!
!    nradfil    number of radar files
!    radfname   file name for radar datasets
!    srcrad     name of radar sources
!
!  OUTPUT:
!
!    isrcrad  index of radar source
!    stnrad   radar site name    character*4
!    latrad   latitude of radar  (degrees N)
!    lonrad   longitude of radar (degrees E)
!    elvrad   elevation of feed horn of radar (m MSL)
!    latradc  latitude of radar column   (degrees N)
!    lonradc  longitude of radar column  (degrees E)
!    irad     radar number of each column
!    nlevrad  number of levels of radar data in each column
!    hgtradc  height (m MSL) of radar observations
!    obsrad   radar observations
!    ncolrad  number of radar columns read-in, local for MPI
!    ncolrad_mpi  number of radar columns read-in, local for MPI
!    istatus  status indicator
!    istat_radar status indicator
!    refmax_mos_3d 3d gridded mosaic of reflectivity (max)
!    ref_mos_3d 3d gridded mosaic of reflectivity (nearest radar)
!    rhv_mos_3d 3d gridded mosaic of Rho-HV, HV Correlation (nearest radar)
!    zdr_mos_3d 3d gridded mosaic of Zdr (nearest radar)
!    kdp_mos_3d 3d gridded mosaic of Kdp (nearest radar)
!
!    tem1     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'mp.inc'
  INCLUDE 'grid.inc'

  INTEGER :: nx,ny,nz
  INTEGER :: nsrc_rad,nvar_radin,mx_rad,nz_rdr,mx_colrad
  INTEGER :: mx_pass
!
  INTEGER :: raduvobs
  INTEGER :: radrhobs
  INTEGER :: radistride
  INTEGER :: radkstride
  INTEGER :: cloudopt
  INTEGER :: iuserad(0:nsrc_rad,mx_pass)
  INTEGER :: iusechk
  INTEGER :: npass
  INTEGER :: np
  REAL    :: xs(nx)
  REAL    :: ys(ny)
  REAL    :: zs(nx,ny,nz)
  REAL    :: latgr(nx,ny)
  REAL    :: longr(nx,ny)
!
  INTEGER :: nradfil
  CHARACTER (LEN=256) :: radfname(mx_rad)
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: srcrad(nsrc_rad)
  INTEGER :: isrcrad(mx_rad)
  REAL :: latrad(mx_rad),lonrad(mx_rad)
  REAL :: elvrad(mx_rad)
  CHARACTER (LEN=5) :: stnrad(mx_rad)
!
!-----------------------------------------------------------------------
!
!  Radar observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: irad(mx_colrad)
  INTEGER :: nlevrad(mx_colrad)
  REAL    :: latradc(mx_colrad),lonradc(mx_colrad)
  REAL    :: hgtradc(nz_rdr,mx_colrad)
  REAL    :: obsrad(nvar_radin,nz_rdr,mx_colrad)
  INTEGER :: indexrad(mx_colrad)
  INTEGER :: oindexrad(mx_colrad)
  REAL    :: refmax_mos_3d(nx,ny,nz)
  REAL    :: ref_mos_3d(nx,ny,nz)
  REAL    :: rhv_mos_3d(nx,ny,nz)
  REAL    :: zdr_mos_3d(nx,ny,nz)
  REAL    :: kdp_mos_3d(nx,ny,nz)
  REAL    :: dist_mos_3d(nx,ny,nz)
  INTEGER :: ncolrad, ncolrad_mpi
  INTEGER :: istatus
  INTEGER :: istat_radar
  INTEGER :: iprocv,jprocv
  INTEGER :: iradvel
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
  REAL :: tem3(nx,ny,nz)
  REAL :: tem4(nx,ny,nz)
  REAL :: tem5(nx,ny,nz)
  REAL :: tem6(nx,ny,nz)
  REAL :: tem7(nx,ny,nz)
  REAL :: tem8(nx,ny,nz)
  REAL :: tem9(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  hdf variables and temporary arrays
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: mxradvr=10
  INTEGER :: iradvr(mxradvr)
!
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp1(:,:,:)
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp2(:,:)
  REAL, ALLOCATABLE :: hmax(:)
  REAL, ALLOCATABLE :: hmin(:)
  REAL, ALLOCATABLE :: latradt(:,:)
  REAL, ALLOCATABLE :: lonradt(:,:)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: velchek=-200.
  REAL, PARAMETER :: refchek=-20.
!
  INTEGER, PARAMETER :: nsrcradin=5
  CHARACTER (LEN=8) :: srcradin(nsrcradin)
  DATA srcradin /'88D-AII','88D-NIDS','CASA-IP1',     &
                 'TDWR','88D-POL'/
!
  CHARACTER (LEN=4)   :: stn
  CHARACTER (LEN=80)  :: runname
  CHARACTER (LEN=256) :: fname
  INTEGER :: ireftim,itime,vcpnum,isource,idummy
  INTEGER :: iradfmt,strhoptin,mapprin,dmpfmt
  INTEGER :: nchanl,ierr,nradvr,ipass, ncolbgn
  INTEGER :: iyr, imon, idy, ihr, imin, isec
  INTEGER :: i,j,k,ivar,isrc,icstrt,icol,jcol
  INTEGER :: klev,kk,kradf,nfile,maxk,knt,lens,sd_id
  INTEGER :: istart,iend,jstart,jend
  REAL    :: xradc,yradc

  REAL :: dxin,dyin,dzin,dzminin,ctrlatin
  REAL :: ctrlonin,tlat1in,tlat2in,tlonin,scalin,rdummy
  REAL :: latcin,loncin
  REAL :: head,sfcrng,elev,srange
  REAL :: xrd,yrd,dummy

  LOGICAL :: hdf_alloc, fndsrcrad, proccols

  INTEGER :: num_colrad, num_verlev
  INTEGER :: irange,iproc,jproc
  INTEGER :: typelev
  REAL    :: xmin, xmax, ymin, ymax
  INTEGER :: numradcol, nummaxlev, nlevsrd

  INTEGER, ALLOCATABLE :: coli(:)
  INTEGER, ALLOCATABLE :: colj(:)
  INTEGER, ALLOCATABLE :: colk(:,:)
  INTEGER, ALLOCATABLE :: numlev(:)
  REAL,    ALLOCATABLE :: collat(:)
  REAL,    ALLOCATABLE :: collon(:)
  REAL,    ALLOCATABLE :: radcolhgt(:,:)
  REAL,    ALLOCATABLE :: radcolref(:,:)
  REAL,    ALLOCATABLE :: radcolrhv(:,:)
  REAL,    ALLOCATABLE :: radcolzdr(:,:)
  REAL,    ALLOCATABLE :: radcolkdp(:,:)
  REAL,    ALLOCATABLE :: radcolvel(:,:)
  REAL,    ALLOCATABLE :: radcolnyq(:,:)
  REAL,    ALLOCATABLE :: radcoltim(:,:)

  INTEGER :: nxlg, nylg
  INTEGER :: ii, jj, isub0, isub, jsub0, jsub
  INTEGER :: kcol
  INTEGER :: tcol
  INTEGER :: nold
  INTEGER :: nnew
  INTEGER :: imax

  REAL    :: refelvmin,refelvmax
  INTEGER :: irngmin,irngmax
  REAL    :: rngmin,rngmax

  REAL, PARAMETER :: epsdx = 0.1
  REAL, PARAMETER :: epsdz = 0.01
  REAL, PARAMETER :: epslat = 0.001

  INTEGER :: dualpol
  LOGICAL :: matched,dualdata
  LOGICAL :: verbose = .TRUE.

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  maxk=0
  icol=0
  tcol=0
  ncolrad_mpi = 0
  istatus=0
  istat_radar=0
  icstrt=0
  nfile=nradfil
  hdf_alloc=.FALSE.
  nold = 0
  nnew = 0
  proccols=(raduvobs > 0 .OR. radrhobs > 0)

  IF(nradfil == 0 ) THEN
    ncolrad = 0
  ELSE IF(nradfil > mx_rad) THEN
    IF (myproc == 0)                                                    &
      WRITE(6,'(a,i3,a,i3/a,i3,a)')                                     &
        ' WARNING nradfil ',nradfil,' exceeds mx_rad dimension ',       &
        mx_rad,' only ',mx_rad,' files will be read.'
    nfile=mx_rad
  END IF

  IF (mp_opt > 0) THEN
    istart = (loc_x-iprocv-1) * (nx-3) + 1  ! this is the global index
    IF (istart < 1 ) istart = 1             ! so the radar grid and the analysis grid
    iend = (loc_x+iprocv) * (nx-3) + 3      ! must be the same, and it also cannot
    k = (nx-3) * nproc_x + 3                ! read split files.
    IF (iend > k) iend = k

    jstart = (loc_y-jprocv-1) * (ny-3) + 1
    IF (jstart < 1 ) jstart = 1
    jend = (loc_y+jprocv) * (ny-3) + 3
    k = (ny-3) * nproc_y + 3
    IF (jend > k) jend = k

  ELSE
    istart = 1
    iend   = nx
    jstart = 1
    jend   = ny
  END IF

  dist_mos_3d=41.E06      ! circumference of earth (m)

!
!-----------------------------------------------------------------------
!
!  Loop through all radars
!
!-----------------------------------------------------------------------
!
  DO kradf=1,nfile
    fname=radfname(kradf)
    dualdata=.FALSE.
    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(fname, '-F f77 -N ieee', ierr)


    lens=LEN(trim(fname))
    IF (fname(lens-3:lens) == 'hdf4') THEN
      dmpfmt=3
    ELSE
      dmpfmt=1
    ENDIF

    IF (myproc == 0) THEN
      WRITE(6,'(/1x,a,a)') 'RDRADCOL: Opening radar file - ',trim(fname)

      IF (dmpfmt == 1) THEN

        CALL getunit( nchanl )
        OPEN(UNIT=nchanl,FILE=trim(fname),ERR=400,                      &
                         FORM='unformatted',STATUS='old')

        istatus=0

        READ(nchanl) stn
        READ(nchanl) ireftim, itime,vcpnum,isource,dualpol,             &
                      idummy,idummy,idummy, idummy,idummy
      ELSE

        CALL hdfopen(trim(fname), 1, sd_id)
        IF (sd_id < 0) THEN
          WRITE (6,'(1x,3a)') 'RDRADCOL: ERROR opening hdf file:',      &
                    trim(fname),' for reading.'
          istatus = sd_id
          GOTO 700
        END IF
        CALL hdfrdc(sd_id, 4, 'radid', stn, istatus)

!
!  Check for corrupt file.
!
        IF (istatus /= 0 ) THEN
          WRITE(6,'(1x,a,1x,a)') 'RDRADCOL:  ERROR corrupt file: ',    &
                    trim(fname),' for reading.'
          CALL hdfclose(sd_id,istatus)
          GOTO 700
        END IF

        CALL hdfrdi(sd_id, 'ireftim', ireftim, istatus)
        CALL hdfrdi(sd_id, 'itime',   itime,   istatus)
        CALL hdfrdi(sd_id, 'vcpnum',  vcpnum,  istatus)
        CALL hdfrdi(sd_id, 'isource', isource, istatus)
        CALL hdfrdi(sd_id, 'dualpol', dualpol, istatus)
      END IF
!
!-----------------------------------------------------------------------
!
!    Connect source number to source name.
!
!-----------------------------------------------------------------------
!
      IF(isource < 1) THEN
        WRITE(6,'(1x,a,i5,/10x,a/)')                                    &
          'WARNING: in rdradcol, read isource as ',isource,             &
                   'Setting isource to default, 1'
        isource=1
      END IF

      IF(isource > nsrcradin) THEN
        WRITE(6,'(1x,a,i3,a,i3,a,/10x,a/)')                             &
          'WARNING: in rdradcol, read isource as ',isource,             &
          ' only ',nsrcradin,' sources known.',                         &
          ' If a new radar source has been added, update rdradcol.f90'
        istatus = -999
        GO TO 700
      END IF

      fndsrcrad = .FALSE.
      DO isrc=1,nsrc_rad
        IF(srcrad(isrc) == srcradin(isource)) THEN
          isrcrad(kradf)=isrc
          WRITE(6,'(1x,a,i3,2a)')                                    &
               'RADAR file source, srcrad(',isrc,') = ',srcrad(isrc)
          fndsrcrad = .TRUE.
          EXIT
        END IF
      END DO

      IF (.NOT. fndsrcrad) THEN
        WRITE(6,'(/1x,a,i4,a,a)') 'ERROR: rdradcol read isource: ',     &
               isource, ' could not find srcrad= ',srcradin(isource)
        WRITE(6,'(1x,a//)') 'Please check radar source names'

        istatus = -999
        GO TO 700
      END IF
!
!-----------------------------------------------------------------------
!
!    Read more header data
!
!-----------------------------------------------------------------------
!
      IF (dmpfmt == 1) THEN

        READ(nchanl) runname
        READ(nchanl) iradfmt,strhoptin,mapprin,irngmin,irngmax,         &
                     typelev,numradcol,nummaxlev,idummy,idummy
        READ(nchanl) dxin,dyin,dzin,dzminin,ctrlatin,                   &
                     ctrlonin,tlat1in,tlat2in,tlonin,scalin,            &
                     latrad(kradf),lonrad(kradf),elvrad(kradf),         &
                     refelvmin,refelvmax
        READ(nchanl) nradvr,iradvr

      ELSE

        CALL hdfrdc(sd_id, 40, 'runname', runname, istatus)
        CALL hdfrdi(sd_id, 'iradfmt', iradfmt, istatus)
        CALL hdfrdi(sd_id, 'strhopt', strhoptin, istatus)
        CALL hdfrdi(sd_id, 'mapproj', mapprin, istatus)
        CALL hdfrdi(sd_id, 'irngmin', irngmin, istatus)
        CALL hdfrdi(sd_id, 'irngmax', irngmax, istatus)

        CALL hdfrdr(sd_id, 'dx', dxin, istatus)
        CALL hdfrdr(sd_id, 'dy', dyin, istatus)
        CALL hdfrdr(sd_id, 'dz', dzin, istatus)
        CALL hdfrdr(sd_id, 'dzmin', dzminin, istatus)
        CALL hdfrdr(sd_id, 'ctrlat', ctrlatin, istatus)
        CALL hdfrdr(sd_id, 'ctrlon', ctrlonin, istatus)
        CALL hdfrdr(sd_id, 'trulat1', tlat1in, istatus)
        CALL hdfrdr(sd_id, 'trulat2', tlat2in, istatus)
        CALL hdfrdr(sd_id, 'trulon', tlonin, istatus)
        CALL hdfrdr(sd_id, 'sclfct', scalin, istatus)
        CALL hdfrdr(sd_id, 'latrad', latrad(kradf), istatus)
        CALL hdfrdr(sd_id, 'lonrad', lonrad(kradf), istatus)
        CALL hdfrdr(sd_id, 'elvrad', elvrad(kradf), istatus)
        CALL hdfrdr(sd_id, 'refelvmin', refelvmin, istatus)
        CALL hdfrdr(sd_id, 'refelvmax', refelvmax, istatus)

        CALL hdfrdi(sd_id, 'nradvr', nradvr, istatus)
        CALL hdfrd1di(sd_id,'iradvr', mxradvr,iradvr,istatus)

        CALL hdfrdi(sd_id, 'typelev', typelev, istatus)
        IF (istatus /= 0) THEN
          typelev = 0
          istatus = 0
        END IF

      END IF

      IF (verbose ) THEN
        WRITE(6,'(1x,2a)')      'runname: ',runname
        WRITE(6,'(1x,3(a,I2))') 'iradfmt: ',iradfmt,                      &
                      ', strhoptin: ',strhoptin,', mapprin: ',mapprin
        WRITE(6,'(1x,a,3F8.0)') 'dxin,dyin,dzin: ',dxin,dyin,dzin
        WRITE(6,'(1x,a,2F8.2)') 'ctrlatin,ctrlonin: ',ctrlatin,ctrlonin
        WRITE(6,'(1x,a,3F8.2)') 'tlat1in,tlat2in,tlonin: ',               &
                                 tlat1in,tlat2in,tlonin
        WRITE(6,'(1x,a,F8.0)')  'scalin: ',scalin
        WRITE(6,'(1x,a,3F8.2)') 'latrad,lonrad,elvrad: ',                 &
                                 latrad(kradf),lonrad(kradf),elvrad(kradf)
        WRITE(6,'(1x,a,20I4)')  'Got nradvr,iradvr: ',nradvr,iradvr
        WRITE(6,'(1x,a,i6)') ' Got dualpol: ',dualpol

        CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
        iyr=MOD(iyr,100)
        WRITE(6,'(/a,i2.2,a,i2.2,a,i2.2,1X,i2.2,a,i2.2,a)')               &
            ' Reading remapped raw radar data for: ',                     &
            imon,'-',idy,'-',iyr,ihr,':',imin,' UTC'

      END IF

      700 CONTINUE

    END IF

    CALL mpbarrier
    CALL mpupdatei(istatus,1)
    IF (istatus == -999) CALL arpsstop(' ',1)
    IF (istatus < 0) CYCLE

    CALL mpupdatec(stn,4)
    stnrad(kradf) = stn

    CALL mpupdatei(dualpol,1)
    CALL mpupdatei(mapprin,1)
    CALL mpupdatei(strhoptin,1)
    CALL mpupdater(dxin,1)
    CALL mpupdater(dyin,1)
    CALL mpupdater(dzin,1)
    CALL mpupdater(dzminin,1)

    CALL mpupdater(ctrlatin,1)
    CALL mpupdater(ctrlonin,1)

    CALL mpupdater(tlat1in,1)
    CALL mpupdater(tlat2in,1)
    CALL mpupdater(tlonin,1)

    CALL mpupdatei(isource,1)
    CALL mpupdatei(typelev,1)
    CALL mpupdatei(irngmin,1)
    CALL mpupdatei(irngmax,1)
    CALL mpupdater(refelvmin,1)
    CALL mpupdater(refelvmax,1)

    CALL mpupdatei(isrcrad(kradf),1)
    CALL mpupdater(latrad (kradf),1)
    CALL mpupdater(lonrad (kradf),1)
    CALL mpupdater(elvrad (kradf),1)
!
!-----------------------------------------------------------------------
!
!   If this radar source is not going to be used in the successive
!   corrections analysis step and cloud analysis processing is off,
!   then don't bother reading columns.
!
!-----------------------------------------------------------------------
!
    dualdata=(dualpol > 0)
    IF(iusechk > 0 ) THEN
      knt=0
      DO ipass=1,npass
        knt=knt+iuserad(isource,ipass)
      END DO
    ELSE
      knt=1
    END IF

    IF (typelev == 0) nold = nold + 1
    IF (typelev == 1) nnew = nnew + 1

    IF (nold > 0 .AND. nnew > 0) THEN
      IF (myproc == 0) WRITE(6,'(1x,a)') 'Mixing old and new radar files is not allowed.'
      CALL arpsstop('Mixing old and new radar files is not allowed.',1)
    END IF

    IF (mp_opt > 0 .AND. nold > 0) THEN
      IF (myproc == 0) WRITE(6,'(1x,a)') 'MPI using radial velocity requires new format data'
      CALL arpsstop('MPI using radial velocity requires new format data',1)
    END IF

    IF( knt > 0 .OR. cloudopt > 0 ) THEN

      IF (typelev == 1) THEN  ! new format with ARPS vertical heights

        IF (myproc == 0) THEN
          IF (dmpfmt == 1) THEN
            !READ(nchanl) xmin, xmax, ymin, ymax  ! Binary does not write those variables
          ELSE
            CALL hdfrdi(sd_id, 'numradcol', numradcol, istatus)
            CALL hdfrdi(sd_id, 'nummaxelv', nummaxlev, istatus)

            CALL hdfrdr(sd_id, 'xmin',   xmin,    istatus)
            CALL hdfrdr(sd_id, 'xmax',   xmax,    istatus)
            CALL hdfrdr(sd_id, 'ymin',   xmin,    istatus)
            CALL hdfrdr(sd_id, 'ymax',   ymax,    istatus)
          END IF
        END IF
        CALL mpupdatei(numradcol,1)
        CALL mpupdatei(nummaxlev,1)

        IF (proccols .AND. nummaxlev > nz_rdr) THEN
          IF (myproc == 0) WRITE(6,'(1x,a,I4,a,/,a,/)')                 &
            'Radar column maximum vertical level nz_rdr (',nz_rdr,') is too small.',  &
            'Please change the value in adas.inc and recompile the program.'
          CALL mpbarrier
          CALL arpsstop('ERROR: Dimension nz_rdr in adas.inc is too small.',1)
        END IF

        ALLOCATE(coli(numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:coli")

        ALLOCATE(colj(numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:colj")

        ALLOCATE(colk(nummaxlev,numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:colk")

        ALLOCATE(numlev(numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:numlev")

        coli(:) = 0
        colj(:) = 0
        colk(:,:) = 0
        numlev(:) = 0

        ALLOCATE(collat(numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:collat")

        ALLOCATE(collon(numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:collon")

        collat(:) = 0.
        collon(:) = 0.

        ALLOCATE(radcolhgt(nummaxlev,numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:radcolhgt")

        ALLOCATE(radcolref(nummaxlev,numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:radcolref")

        ALLOCATE(radcolvel(nummaxlev,numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:radcolvel")

        ALLOCATE(radcolnyq(nummaxlev,numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:radcolnyq")

        ALLOCATE(radcoltim(nummaxlev,numradcol), STAT = istatus)
        CALL check_alloc_status(istatus, "rdradcol:radcoltim")

        radcolhgt(:,:) = -999.
        radcolref(:,:) = -999.
        radcolvel(:,:) = -999.
        radcolnyq(:,:) = -999.
        radcoltim(:,:) = -999.

        IF(dualdata) THEN
          ALLOCATE(radcolrhv(nummaxlev,numradcol), STAT = istatus)
          CALL check_alloc_status(istatus, "rdradcol:radcolrhv")

          ALLOCATE(radcolzdr(nummaxlev,numradcol), STAT = istatus)
          CALL check_alloc_status(istatus, "rdradcol:radcolzdr")

          ALLOCATE(radcolkdp(nummaxlev,numradcol), STAT = istatus)
          CALL check_alloc_status(istatus, "rdradcol:radcolkdp")
          radcolrhv(:,:) = -999.
          radcolzdr(:,:) = -999.
          radcolkdp(:,:) = -999.
        END IF

        IF (myproc == 0) THEN
          IF (dmpfmt == 1) THEN ! READ binary files

            kcol  = 0
            DO jcol=1,numradcol
              READ(nchanl,END=201) i,j,xrd,yrd,                         &
                                   latcin,loncin,elev,klev

              kcol = kcol + 1
              coli(kcol) = i
              colj(kcol) = j

              collat(kcol)  = latcin
              collon(kcol)  = loncin
              numlev(kcol)  = klev

              READ(nchanl,END=202)      (colk(kk,kcol),kk=1,klev)
              READ(nchanl,END=202) (radcolhgt(kk,kcol),kk=1,klev)
              READ(nchanl,END=202) (radcolref(kk,kcol),kk=1,klev)
              READ(nchanl,END=202) (radcolvel(kk,kcol),kk=1,klev)
              READ(nchanl,END=202) (radcolnyq(kk,kcol),kk=1,klev)
              READ(nchanl,END=202) (radcoltim(kk,kcol),kk=1,klev)
              IF(dualdata) THEN
                READ(nchanl,END=201) (radcolrhv(kk,kcol),kk=1,klev)
                READ(nchanl,END=201) (radcolzdr(kk,kcol),kk=1,klev)
                READ(nchanl,END=201) (radcolkdp(kk,kcol),kk=1,klev)
              END IF
            END DO

            201   CONTINUE
            WRITE(6,'(a,i6,a)') ' End of file reached after reading',     &
                             (kcol-icstrt),' columns'
            GO TO 205

            202   CONTINUE
            WRITE(6,'(a,i6,a)') ' End of file reached while reading',     &
                                kcol,' column'

            205   CONTINUE

            CLOSE(nchanl)
            CALL retunit( nchanl )

            nummaxlev = MAX(nummaxlev,klev)

          ELSE                  ! READ hdf file

            ALLOCATE (itmp2(nummaxlev,numradcol),stat=istatus)
            CALL check_alloc_status(istatus, "rdradcol:itmp2")

            CALL hdfrd1d (sd_id,'radcollat',numradcol,collat,istatus)
            CALL hdfrd1d (sd_id,'radcollon',numradcol,collon,istatus)
            CALL hdfrd1di(sd_id,'numelev',  numradcol,numlev,istatus)

            CALL hdfrd1di(sd_id,'radcoli',  numradcol,coli,istatus)
            CALL hdfrd1di(sd_id,'radcolj',  numradcol,colj,istatus)

            CALL hdfrd2di(sd_id,'radcolk',  nummaxlev,numradcol,colk,istatus)

            CALL hdfrd2d (sd_id,'radcolhgt',nummaxlev,numradcol,radcolhgt,istatus,itmp2)
            CALL hdfrd2d (sd_id,'radcolref',nummaxlev,numradcol,radcolref,istatus,itmp2)
            CALL hdfrd2d (sd_id,'radcolvel',nummaxlev,numradcol,radcolvel,istatus,itmp2)
            CALL hdfrd2d (sd_id,'radcolnyq',nummaxlev,numradcol,radcolnyq,istatus,itmp2)
            CALL hdfrd2d (sd_id,'radcoltim',nummaxlev,numradcol,radcoltim,istatus,itmp2)
            IF(dualdata) THEN
              CALL hdfrd2d (sd_id,'radcolrhv',nummaxlev,numradcol,radcolrhv,istatus,itmp2)
              CALL hdfrd2d (sd_id,'radcolzdr',nummaxlev,numradcol,radcolzdr,istatus,itmp2)
              CALL hdfrd2d (sd_id,'radcolkdp',nummaxlev,numradcol,radcolkdp,istatus,itmp2)
            END IF

            CALL hdfclose(sd_id,istatus)

            DEALLOCATE(itmp2)
          END IF

        END IF    ! myproc == 0

        CALL mpupdatei(coli,numradcol)
        CALL mpupdatei(colj,numradcol)
        CALL mpupdatei(colk,numradcol*nummaxlev)
        CALL mpupdatei(numlev,numradcol)
        CALL mpupdater(collat,numradcol)
        CALL mpupdater(collon,numradcol)
        CALL mpupdater(radcolhgt,nummaxlev*numradcol)
        CALL mpupdater(radcolref,nummaxlev*numradcol)
        CALL mpupdater(radcolvel,nummaxlev*numradcol)
        CALL mpupdater(radcolnyq,nummaxlev*numradcol)
        CALL mpupdater(radcoltim,nummaxlev*numradcol)
        IF (dualdata) THEN
          CALL mpupdater(radcolzdr,nummaxlev*numradcol)
          CALL mpupdater(radcolrhv,nummaxlev*numradcol)
          CALL mpupdater(radcolkdp,nummaxlev*numradcol)
        END IF

        IF(proccols) THEN
          DO kcol = 1, numradcol
            IF (MOD(coli(kcol)+colj(kcol),radistride) == 0) THEN
              tcol = tcol + 1           ! absolute column number
              IF (iradvel > 0) THEN
                IF (coli(kcol) < istart .OR. coli(kcol) > iend) CYCLE
                IF (colj(kcol) < jstart .OR. colj(kcol) > jend) CYCLE
              END IF
              icol = icol + 1
              IF (icol > mx_colrad ) THEN
                 WRITE(6,'(1x,a,I5,a,/,a,/)')   &
             'Radar column number mx_colrad (',mx_colrad,') is too small.',  &
             'Please change the value in adas.inc and recompile the program.'
                 CALL arpsstop( &
                      'ERROR: Dimension mx_colrad in adas.inc is too small.',1)
              END IF

              IF (mp_opt > 0) THEN
                CALL lltoxy(1,1,collat(kcol),collon(kcol),xradc,yradc)
                indexrad(icol)  = -1
                oindexrad(icol) = tcol
                IF( xradc >= xs(1) .AND. xradc <= xs(nx-1)  .AND.        &
                    yradc >= ys(1) .AND. yradc <= ys(ny-1)) THEN
                     indexrad(icol) = myproc   ! temporary setting for using
                                               ! with in the call of prepradar
                ENDIF
              ENDIF
              irad(icol)    = kradf
              latradc(icol) = collat(kcol)
              lonradc(icol) = collon(kcol)
              nlevrad(icol) = numlev(kcol)
              DO kk = 1,numlev(kcol)
                IF ( MOD( (colk(kk,kcol)-1), radkstride ) == 0  .OR.      &
                     kk == numlev(kcol) ) THEN
                  hgtradc (kk,icol) = radcolhgt(kk,kcol)
                  obsrad(1,kk,icol) = radcolref(kk,kcol)
                  obsrad(2,kk,icol) = radcolvel(kk,kcol)
                  obsrad(3,kk,icol) = radcolnyq(kk,kcol)
                  obsrad(4,kk,icol) = radcoltim(kk,kcol)
                END IF
              END DO
            END IF  ! horizontal stride
          END DO

          IF ( verbose .AND. (icol-icstrt) > 0 ) THEN
            WRITE(6,'(1x,a,I5,a,i6,3a)') 'Processor: ',myproc,' accepts ',&
            (icol-icstrt),' non-missing columns from file - ',trim(fname),'.'
          END IF

          icstrt = icol
          maxk   = nummaxlev

        END IF   !  proccols
!
! Process for radar mosaic
!
        IF(cloudopt > 0) THEN
!
!  First check for match of remapped grid and analysis grid
!
          matched=.TRUE.
          IF( mapprin /= mapproj .OR. strhoptin /= strhopt ) THEN
            matched=.FALSE.
            WRITE(6,'(1x,I3,a,2a,/2(a,i10))') myproc,':',               &
              'WARNING: Grid mis-match for remapped radar ',stn,        &
              '  Radar mapproj:',mapprin,'  strhopt:',strhoptin
          END IF
          IF(abs(dxin-dx) > epsdx .OR. abs(dyin-dy) > epsdx .OR.        &
             abs(dzin-dz) > epsdz .OR. abs(dzminin-dzmin) > epsdz ) THEN
            matched=.FALSE.
            WRITE(6,'(1x,I3,3a,/2(a,f12.1),2(a,f12.3))') myproc,':',    &
              'WARNING: Grid mis-match for remapped radar ',stn,        &
              '  Radar dx:',dxin,'  dy:',dyin, &
              '    dz:',dzin,'  dzmin:',dzminin
          END IF
          IF(abs(ctrlatin-ctrlat) > epslat .OR.                         &
             abs(ctrlonin-ctrlon) > epslat ) THEN
            matched=.FALSE.
            WRITE(6,'(1x,I3,3a,/2(a,f12.4))') myproc, ':',              &
              'WARNING: Grid mis-match for remapped radar ',stn,        &
              '  Radar ctrlat:',ctrlatin,'  ctrlonin:',ctrlonin
          END IF
          IF(abs(tlat1in-trulat1) > epslat .OR.                         &
             abs(tlat2in-trulat2) > epslat .OR.                         &
             abs(tlonin  -trulon) > epslat ) THEN
            matched=.FALSE.
            WRITE(6,'(1x,I3,3a,/3(a,f12.4))') myproc, ':',              &
              'WARNING: Grid mis-match for remapped radar ',stn,        &
              '  Radar trulat1:',tlat1in,'  trulat2:',tlat2in,          &
              '    trulon:',tlonin
          END IF
!
!  Since clear columns are not written, set reflectivity mosaic
!  to zero within the active radar scan area.  The zero value replaces
!  the default missing value and carries a different meaning in the
!  cloud analysis.  Like the mosaicking done for observed reflectivity
!  values this is done with a max function to avoid overwriting data
!  from previous radars.
!
          IF(matched) THEN
            print *, ' MATCHED !!! '
            IF(myproc == 0) THEN
              WRITE(6,'(1x,3a)') 'Adding radar ',stn,' to mosaic.'
              WRITE(6,'(1x,2(a,f10.2),a)') 'Reflectivity coverage ',    &
                (irngmin*0.001),' km  to ',(irngmax*0.001),' km.'
              WRITE(6,'(1x,2(a,f10.2),a)') 'Reflectivity elevs    ',    &
                refelvmin,' deg to ',refelvmax,' deg.'
            END IF
            DO j=1,ny
              DO i=1,nx
                CALL disthead(latrad(kradf),lonrad(kradf),              &
                              latgr(i,j),longr(i,j),head,sfcrng)
                DO k=2,nz-1
                  CALL beamelv((zs(i,j,k)-elvrad(kradf)),sfcrng,elev,srange)
                  irange = NINT(srange)
                  IF(irange < irngmax  .AND. irange > irngmin .AND.     &
                     elev >= refelvmin .AND. elev <= refelvmax ) THEN
                    refmax_mos_3d(i,j,k)=max(refmax_mos_3d(i,j,k),0.)
                    ref_mos_3d(i,j,k)=max(ref_mos_3d(i,j,k),0.)
                  END IF
                END DO
              END DO
            END DO
!
!  Now process the radar columns that are in this processor's domain.
!  ref_max mosaic is formed using the max value of all radar observations.
!  Mosaic of ref and dual pol variables is formed using the nearest radar observations.
!
            IF(dualdata) THEN
              DO kcol = 1, numradcol
                iproc = (coli(kcol)-2)/(nx-3) + 1
                jproc = (colj(kcol)-2)/(ny-3) + 1

                IF (iproc == loc_x .AND. jproc == loc_y) THEN
                  i =  MOD((coli(kcol)-2),(nx-3)) + 2
                  j =  MOD((colj(kcol)-2),(ny-3)) + 2

                  CALL disthead(latrad(kradf),lonrad(kradf),              &
                                latgr(i,j),longr(i,j),head,sfcrng)

                  DO kk=1,numlev(kcol)
                    k=colk(kk,kcol)
                    IF(k <= nz .AND. k > 0) THEN
                      refmax_mos_3d(i,j,k)=max(refmax_mos_3d(i,j,k),radcolref(kk,kcol))
                      IF(sfcrng < dist_mos_3d(i,j,k) .and. radcolref(kk,kcol) > refchek) THEN
                        dist_mos_3d(i,j,k)=sfcrng
                        ref_mos_3d(i,j,k)=radcolref(kk,kcol)
                        rhv_mos_3d(i,j,k)=radcolrhv(kk,kcol)
                        zdr_mos_3d(i,j,k)=radcolzdr(kk,kcol)
                        kdp_mos_3d(i,j,k)=radcolkdp(kk,kcol)
                      END IF
                    END IF  ! 1 < k < nz
                  END DO  ! kk = 1, klev
                END IF  ! this processor
              END DO  ! kcol
            ELSE  ! no dual-pol data
              DO kcol = 1, numradcol
                iproc = (coli(kcol)-2)/(nx-3) + 1
                jproc = (colj(kcol)-2)/(ny-3) + 1

                IF (iproc == loc_x .AND. jproc == loc_y) THEN
                  i =  MOD((coli(kcol)-2),(nx-3)) + 2
                  j =  MOD((colj(kcol)-2),(ny-3)) + 2

                  IF( i < (nx-2) .AND. j < (nx-2) .AND.                &
                      i > 0 .AND. j > 0 ) THEN

!                   print *, ' coli,colj: ',coli(kcol),colj(kcol)

                    CALL disthead(latrad(kradf),lonrad(kradf),           &
                                  latgr(i,j),longr(i,j),head,sfcrng)

                    DO kk=1,numlev(kcol)
                      k=colk(kk,kcol)
                      IF(k <= nz .AND. k > 0) THEN
                        refmax_mos_3d(i,j,k)=max(refmax_mos_3d(i,j,k),radcolref(kk,kcol))
                        IF(sfcrng < dist_mos_3d(i,j,k) .and. radcolref(kk,kcol) > refchek) THEN
                          dist_mos_3d(i,j,k)=sfcrng
                          ref_mos_3d(i,j,k)=radcolref(kk,kcol)
                        END IF
                      END IF  ! 1 < k < nz
                    END DO  ! kk = 1, klev
                  END IF  ! col index within limits
                END IF  ! this processor
              END DO  ! kcol
            END IF  ! dualdata
            istat_radar=1
          END IF ! grid matched

        END IF ! cloudopt ref_mos_3d_block

        DEALLOCATE(coli, colj, colk,numlev)
        DEALLOCATE(collat, collon)

        DEALLOCATE(radcolhgt, radcolref, radcolvel, radcolnyq, radcoltim)

      ELSE       ! Maybe old format
!
!-----------------------------------------------------------------------
!
!  Note here the radar data indices:
!
!       1 Reflectivity
!       2 Radial Velocity
!       3 Nyquist Velocity
!       4 Time in seconds from the reference time
!
!-----------------------------------------------------------------------
!
        IF (dmpfmt == 1) THEN
          DO jcol=1,999999
            READ(nchanl,END=101) i,j,xrd,yrd,                           &
                latcin,loncin,elev,klev
            READ(nchanl,END=102) (tem1(1,1,kk),kk=1,klev)
            READ(nchanl,END=102) (tem2(1,1,kk),kk=1,klev)
            READ(nchanl,END=102) (tem3(1,1,kk),kk=1,klev)
            READ(nchanl,END=102) (tem4(1,1,kk),kk=1,klev)
            READ(nchanl,END=102) (tem5(1,1,kk),kk=1,klev)
            READ(nchanl,END=102) (tem6(1,1,kk),kk=1,klev)
            READ(nchanl,END=102) (tem7(1,1,kk),kk=1,klev)
            READ(nchanl,END=102) (tem8(1,1,kk),kk=1,klev)
            IF(proccols .AND. mod((i+j),radistride) == 0 ) THEN
              icol=icol+1
              IF(icol <= mx_colrad) THEN
                irad(icol)=kradf
                latradc(icol)=latcin
                lonradc(icol)=loncin
                DO kk=1,nz_rdr
                  hgtradc(kk,icol)=-999.
                  DO ivar=1,nvar_radin
                    obsrad(ivar,kk,icol)=-999.
                  END DO
                END DO
                IF(klev > 1) THEN
                  DO kk=1,(klev-1),radkstride
                    k=nint(tem1(1,1,kk))
                    maxk=MAX(maxk,k)
                    k=MIN(k,nz_rdr)
                    hgtradc(k,icol)=tem2(1,1,kk)
                    obsrad(1,k,icol)=tem3(1,1,kk)
                    obsrad(2,k,icol)=tem4(1,1,kk)
                    obsrad(3,k,icol)=tem5(1,1,kk)
                    obsrad(4,k,icol)=tem6(1,1,kk)
                    nlevrad(icol)=k
                  END DO
                END IF
                k=nint(tem1(1,1,klev))
                maxk=MAX(maxk,k)
                k=MIN(k,nz_rdr)
                hgtradc(k,icol)=tem2(1,1,klev)
                obsrad(1,k,icol)=tem3(1,1,klev)
                obsrad(2,k,icol)=tem4(1,1,klev)
                obsrad(3,k,icol)=tem5(1,1,klev)
                obsrad(4,k,icol)=tem6(1,1,klev)
                nlevrad(icol)=k
              ELSE
                maxk=MAX(maxk,nint(tem1(1,1,klev)))
              END IF
            END IF
!
!  Process for radar mosaic
!
            IF( cloudopt > 0 ) THEN
!
!  Mosaic is formed using the max value of all radar observations.
!  Note: Old format does not include necessary scan info for setting clear.
!
              DO kk=1,klev
                k=nint(tem1(1,1,kk))
                k=MIN(k,nz_rdr)
                ref_mos_3d(i,j,k)=max(ref_mos_3d(i,j,k),tem3(1,1,k))
              END DO

            END IF  ! cloudopt ref_mos_3d block
          END DO
          !istat_radar=1
          101   CONTINUE
          WRITE(6,'(a,i6,a)') ' End of file reached after reading',     &
                             (icol-icstrt),' columns'
          GO TO 105
          102   CONTINUE
          WRITE(6,'(a,i6,a)') ' End of file reached while reading',     &
                             icol,' column'
          icol=icol-1
          105   CONTINUE
          CLOSE(nchanl)
          CALL retunit( nchanl )
          icstrt=icol+1

          IF (icol > 0) istat_radar=1

        ELSE   ! hdf file

          !
          !  Allocate hdf temporary arrays if not already done
          !
          IF( .NOT. hdf_alloc ) THEN
            ALLOCATE(latradt(nx,ny),stat=istatus)
            CALL check_alloc_status(istatus, "rdradcol:latradt")
            latradt=-999999.

            ALLOCATE(lonradt(nx,ny),stat=istatus)
            CALL check_alloc_status(istatus, "rdradcol:lonradt")
            lonradt=-999999.

            ALLOCATE (hmax(nz),stat=istatus)
            CALL check_alloc_status(istatus, "rdradcol:hmax")

            ALLOCATE (hmin(nz),stat=istatus)
            CALL check_alloc_status(istatus, "rdradcol:hmin")

            ALLOCATE (itmp1(nx,ny,nz),stat=istatus)
            CALL check_alloc_status(istatus, "rdradcol:itmp1")

            hdf_alloc=.TRUE.
          END IF

          CALL hdfrd2d(sd_id,'grdlatc',nx,ny,latradt,istatus,itmp1)
          IF (istatus /= 0) GO TO 115
          CALL hdfrd2d(sd_id,'grdlonc',nx,ny,lonradt,istatus,itmp1)
          IF (istatus /= 0) GO TO 115
          CALL hdfrd3d(sd_id,'hgtrad',nx,ny,nz,tem2,istatus,itmp1,hmax,hmin)
          IF (istatus /= 0) GO TO 115
          CALL hdfrd3d(sd_id,'gridref',nx,ny,nz,tem3,istatus,itmp1,hmax,hmin)
          IF (istatus /= 0) GO TO 115
          CALL hdfrd3d(sd_id,'gridvel',nx,ny,nz,tem4,istatus,itmp1,hmax,hmin)
          IF (istatus /= 0) GO TO 115
          CALL hdfrd3d(sd_id,'gridnyq',nx,ny,nz,tem5,istatus,itmp1,hmax,hmin)
          IF (istatus /= 0) GO TO 115
          CALL hdfrd3d(sd_id,'gridtim',nx,ny,nz,tem6,istatus,itmp1,hmax,hmin)
          IF (istatus /= 0) GO TO 115

          IF(proccols) THEN
            DO j=1,ny
              DO i=1,nx
                IF( mod((i+j),radistride) == 0 ) THEN
                  !
                  !  Check for non-missing data in this column
                  !
                  knt=0
                  DO k=1,(nz-1),radkstride
                    IF(tem3(i,j,k) > refchek .OR. &
                       tem4(i,j,k) > velchek ) THEN
                      maxk=MAX(maxk,k)
                      knt=knt+1
                    END IF
                  END DO
                  !
                  !  If non-missing data exist increment column counter and
                  !  transfer to column data
                  !
                  IF( knt > 0 ) THEN
                    icol=icol+1
                    IF(icol <= mx_colrad) THEN
                      irad(icol)=kradf
                      latradc(icol)=latradt(i,j)
                      lonradc(icol)=lonradt(i,j)
                      DO kk=1,nz_rdr
                        hgtradc(kk,icol)=-999.
                        DO ivar=1,nvar_radin
                          obsrad(ivar,kk,icol)=-999.
                        END DO
                      END DO
                      DO k=1,(nz-1),radkstride
                        IF(tem3(i,j,k) > refchek .OR. &
                           tem4(i,j,k) > velchek ) THEN
                          hgtradc(k,icol)=tem2(i,j,k)
                          obsrad(1,k,icol)=tem3(i,j,k)
                          obsrad(2,k,icol)=tem4(i,j,k)
                          obsrad(3,k,icol)=tem5(i,j,k)
                          obsrad(4,k,icol)=tem6(i,j,k)
                          nlevrad(icol)=k
                        END IF
                      END DO
                    END IF  ! icol less than mx_colrad
                  END IF  ! non-missing data in column
                END IF ! stride check
              END DO
            END DO

            IF (myproc == 0) WRITE(6,'(a,i6,a)') ' Read ',(icol-icstrt),&
                  ' non-missing columns from hdf file.'
            icstrt=icol+1
          END IF

!-----------------------------------------------------------------------
!
!  Process for radar mosaic
!
!-----------------------------------------------------------------------
          IF( cloudopt > 0 ) THEN
            !
            !  Now process the radar columns.
            !  Mosaic is formed using the max value of all radar observations.
            !  Note: Old format does not include necessary scan info for setting clear.
            !
            DO j=1,ny
              DO i=1,nx
                CALL disthead(latrad(kradf),lonrad(kradf),              &
                              latgr(i,j),longr(i,j),head,sfcrng)
                DO k=1,nz
                  refmax_mos_3d(i,j,k)=max(refmax_mos_3d(i,j,k),tem3(i,j,k))
                  IF(sfcrng < dist_mos_3d(i,j,k) .and. radcolref(kk,kcol) > refchek) THEN
                    dist_mos_3d(i,j,k)=sfcrng
                    ref_mos_3d(i,j,k)=radcolref(kk,kcol)
                  END IF
                END DO
              END DO
            END DO

            istat_radar=1

          END IF  ! cloudopt ref_mos_3d block

          CALL hdfclose(sd_id,istatus)

        END IF   ! HDF file, old format

      END IF   ! New or Old format?

    ELSE  ! iuserad check failed
      IF( dmpfmt == 1 ) THEN
        IF (mp_opt == 0)                                                &
          WRITE(6,'(a,a,/a/)')                                          &
            '   Data for source ',srcradin(isource),                    &
            ' are not used in successive correction step, skipping'
        CLOSE(nchanl)
        CALL retunit( nchanl )
      ELSE  ! hdf format file
        CALL hdfclose(sd_id,istatus)
      END IF
    END IF  ! use source ?

  END DO  ! file loop

  IF (ALLOCATED(itmp1)) DEALLOCATE(itmp1)

  IF (ALLOCATED(latradt)) THEN
    DEALLOCATE(latradt, lonradt)
    DEALLOCATE(hmax, hmin)
  END IF

  IF(icol > mx_colrad .AND. myproc == 0) THEN
    WRITE(6,'(//a)')       ' #### WARNING max_colrad EXCEEDED ####'
    WRITE(6,'(a,i6,a)')    ' increase mx_colrad from ',mx_colrad,' columns'
    WRITE(6,'(2x,i6,a//)') icol,' columns are needed to read all files'
  END IF
  ncolrad=MIN(icol,mx_colrad)

  IF (mp_opt > 0) ncolrad_mpi = tcol

  IF ( verbose .AND. maxk > 0) THEN         ! Debugging message only
    WRITE(6,'(1x,2(a,i5),2a/)')   &
       'Processor: ',myproc,' finds maximum number of vert levels ',    &
       maxk,' from file - ',TRIM(fname)
    CALL flush(6)
  END IF

  IF (mp_opt > 0) THEN
    imax = icol
    CALL mpmaxi(imax)
    IF (myproc == 0) WRITE(6,'(a,i7/)') ' Maximum local columns ',imax
  END IF

  IF(maxk > nz_rdr) THEN
     IF (myproc == 0)                                                    &
       WRITE(6,'(a,i5)') '   EXCEEDS nz_rdr, increase nz_rdr: ',nz_rdr
     CALL mpbarrier
     CALL arpsstop("nz_rdr is too small",1)
  END IF

  RETURN

  115 CONTINUE
  IF (myproc == 0)                                                      &
    WRITE(6,'(/a/)') ' Error reading data in RDRADCOL(HDF format).'
  CALL mpbarrier
  CALL arpsstop("Problem reading HDF data",1)

  400 CONTINUE
  IF (myproc == 0)                                                      &
    WRITE(6,'(a,a,/a)') '   Error opening radar file ',trim(fname),     &
                      ' Stopping in rdradcol.'
  CALL arpsstop("Error opening file",1)

END SUBROUTINE rdradcol
